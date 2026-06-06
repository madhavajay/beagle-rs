//! Port of `imp/ImputedVcfWriter.java` — writes the observed + imputed genotypes for one target
//! marker cluster (and its surrounding imputed reference markers) as VCF records, combining each
//! haplotype's HMM state probabilities into per-allele probabilities (with linear interpolation
//! between flanking target clusters) and emitting them via `ImputedRecBuilder`.

use std::io::Write;

use crate::blbutil::FloatList;
use crate::ints::IntList;
use crate::vcf::{Samples, GT};

use super::{ImpData, ImputedRecBuilder, RefHapHash, StateProbs};

/// Port of `imp/ImputedVcfWriter.java`.
pub struct ImputedVcfWriter<'a> {
    imp_data: &'a ImpData,
    targ_samples: Samples,
    targ_cluster: i32,
    ref_start: i32,
    clust_end: i32,
    ref_end: i32,
    indices: IntList,
    hashes: IntList,
    seq_probs: FloatList,
    seq_probs_p1: FloatList,
}

impl<'a> ImputedVcfWriter<'a> {
    /// `new ImputedVcfWriter(ImpData impData, int refStart, int refEnd, int targCluster)`.
    pub fn new(imp_data: &'a ImpData, ref_start: i32, ref_end: i32, targ_cluster: i32) -> Self {
        assert!(ref_start >= 0, "{ref_start}");
        assert!(ref_end <= imp_data.ref_gt().n_markers(), "{ref_end}");
        let targ_samples = imp_data.targ_gt().samples().clone();
        let real_ref_start = if targ_cluster == 0 {
            ref_start
        } else {
            ref_start.max(imp_data.ref_cluster_start(targ_cluster))
        };
        let (clust_end, real_ref_end) = if targ_cluster < imp_data.n_clusters() - 1 {
            let tmp_clust_end = ref_start.max(imp_data.ref_cluster_end(targ_cluster));
            (
                tmp_clust_end.min(ref_end),
                imp_data.ref_cluster_start(targ_cluster + 1).min(ref_end),
            )
        } else {
            (ref_end, ref_end)
        };
        ImputedVcfWriter {
            imp_data,
            targ_samples,
            targ_cluster,
            ref_start: real_ref_start,
            clust_end,
            ref_end: real_ref_end,
            indices: IntList::with_capacity(4),
            hashes: IntList::with_capacity(4),
            seq_probs: FloatList::new(),
            seq_probs_p1: FloatList::new(),
        }
    }

    /// `appendRecords(AtomicReferenceArray<StateProbs> stateProbs, PrintWriter out)`.
    pub fn append_records(&mut self, state_probs: &[Box<dyn StateProbs>], out: &mut dyn Write) {
        if self.ref_start >= self.ref_end {
            return;
        }
        let ref_hap_hash = RefHapHash::new(
            state_probs,
            self.targ_cluster,
            self.imp_data.ref_gt().clone(),
            self.ref_start,
            self.ref_end,
        );
        let mut rec_builders = self.rec_builders();
        let mut a1_probs = self.al_probs();
        let mut a2_probs = self.al_probs();
        let is_imputed = self.is_imputed();
        let n = state_probs.len();
        let n_markers = a1_probs.len();
        let mut h = 0;
        while h < n {
            let is_diploid = self.targ_samples.is_diploid((h >> 1) as i32);
            self.set_al_probs(state_probs[h].as_ref(), &ref_hap_hash, &mut a1_probs);
            self.set_al_probs(state_probs[h + 1].as_ref(), &ref_hap_hash, &mut a2_probs);
            for m in 0..n_markers {
                if !is_imputed[m] {
                    self.set_to_obs_alleles(&mut a1_probs, &mut a2_probs, m, h as i32);
                }
                if is_diploid {
                    rec_builders[m].add_sample_data(&mut a1_probs[m], &mut a2_probs[m]);
                } else {
                    rec_builders[m].add_sample_data_haploid(&mut a1_probs[m]);
                }
                a1_probs[m].iter_mut().for_each(|x| *x = 0.0);
                a2_probs[m].iter_mut().for_each(|x| *x = 0.0);
            }
            h += 2;
        }
        for (m, rb) in rec_builders.iter().enumerate() {
            rb.print_rec(out, is_imputed[m]);
        }
    }

    fn set_al_probs(
        &mut self,
        state_probs: &dyn StateProbs,
        rhh: &RefHapHash,
        al_probs: &mut [Vec<f32>],
    ) {
        self.indices.clear();
        self.hashes.clear();
        self.seq_probs.clear();
        self.seq_probs_p1.clear();
        let mut alleles = vec![0i32; (rhh.end() - rhh.start()) as usize];
        for j in 0..state_probs.n_states(self.targ_cluster) {
            let hap = state_probs.ref_hap(self.targ_cluster, j);
            let val = state_probs.probs(self.targ_cluster, j);
            let val_p1 = state_probs.probs_p1(self.targ_cluster, j);
            let index = rhh.hap2_index(hap);
            let hash = rhh.hash(index);
            let mut i = 0;
            while i < self.hashes.size() && self.hashes.get(i) != hash {
                i += 1;
            }
            if i == self.hashes.size() {
                self.indices.add(index);
                self.hashes.add(hash);
                self.seq_probs.add(val);
                self.seq_probs_p1.add(val_p1);
            } else {
                self.seq_probs.add_to_element(i, val);
                self.seq_probs_p1.add_to_element(i, val_p1);
            }
        }
        self.set_al_probs2(al_probs, rhh, &mut alleles);
    }

    fn set_al_probs2(&self, al_probs: &mut [Vec<f32>], rhh: &RefHapHash, alleles: &mut [i32]) {
        let n_seq = self.seq_probs.size();
        if n_seq == 1 {
            let index = self.indices.get(0);
            rhh.set_alleles(index, alleles);
            for m in self.ref_start..self.ref_end {
                let mm = (m - self.ref_start) as usize;
                al_probs[mm][alleles[mm] as usize] = 1.0;
            }
        } else {
            for j in 0..n_seq {
                let index = self.indices.get(j);
                rhh.set_alleles(index, alleles);
                let prob = self.seq_probs.get(j);
                let prob_p1 = self.seq_probs_p1.get(j);
                for m in self.ref_start..self.clust_end {
                    let mm = (m - self.ref_start) as usize;
                    al_probs[mm][alleles[mm] as usize] += prob;
                }
                for m in self.clust_end..self.ref_end {
                    let wt = self.imp_data.weight(m);
                    let mm = (m - self.ref_start) as usize;
                    let a = alleles[mm] as usize;
                    // Java: float[] += (double); float←(float+double) narrowed
                    al_probs[mm][a] = (al_probs[mm][a] as f64
                        + (wt * prob as f64 + (1.0 - wt) * prob_p1 as f64))
                        as f32;
                }
            }
        }
    }

    fn set_to_obs_alleles(
        &self,
        a1_probs: &mut [Vec<f32>],
        a2_probs: &mut [Vec<f32>],
        m: usize,
        targ_hap: i32,
    ) {
        a1_probs[m].iter_mut().for_each(|x| *x = 0.0);
        a2_probs[m].iter_mut().for_each(|x| *x = 0.0);
        let pre_clust_index = self
            .imp_data
            .marker_indices()
            .marker_to_targ_marker_at(self.ref_start + m as i32);
        let a1 = self.imp_data.targ_gt().allele(pre_clust_index, targ_hap);
        let a2 = self
            .imp_data
            .targ_gt()
            .allele(pre_clust_index, targ_hap + 1);
        a1_probs[m][a1 as usize] = 1.0;
        a2_probs[m][a2 as usize] = 1.0;
    }

    fn rec_builders(&self) -> Vec<ImputedRecBuilder> {
        let ref_gt = self.imp_data.ref_gt();
        let gp = self.imp_data.par().gp();
        let ap = self.imp_data.par().ap();
        let n_input_targ_haps = self.imp_data.n_input_targ_haps();
        (self.ref_start..self.ref_end)
            .map(|m| ImputedRecBuilder::new(ref_gt.marker(m).clone(), n_input_targ_haps, ap, gp))
            .collect()
    }

    fn al_probs(&self) -> Vec<Vec<f32>> {
        let ref_markers = self.imp_data.ref_gt().markers();
        (self.ref_start..self.ref_end)
            .map(|m| vec![0f32; ref_markers.marker(m).n_alleles() as usize])
            .collect()
    }

    fn is_imputed(&self) -> Vec<bool> {
        let wm = self.imp_data.marker_indices();
        (0..(self.ref_end - self.ref_start))
            .map(|j| wm.marker_to_targ_marker_at(self.ref_start + j) == -1)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::imp::ImpLS;
    use crate::main_pkg::Par;
    use crate::vcf::{
        allele_ref_gt_rec_from_parser, BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices,
        MarkerParser, PositionMap, RefGT, RefGTRec, VcfHeader, VcfRecGTParser, Window,
        HEADER_PREFIX,
    };
    use std::rc::Rc;

    fn header(n: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n {
            hdr.push_str(&format!("\tS{s}"));
        }
        VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &vec![true; n],
        )
    }

    fn imp_data() -> ImpData {
        let gt = std::env::temp_dir().join("beagle_rs_impvcf_gt.vcf");
        let rf = std::env::temp_dir().join("beagle_rs_impvcf_ref.vcf");
        std::fs::write(&gt, b"x").unwrap();
        std::fs::write(&rf, b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            format!("ref={}", rf.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
        ]);
        let hr = header(6);
        let ht = header(3);
        let mp = MarkerParser::new(true, true, true, true);
        let ref_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|1\t1|0\t0|0",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\t1|1\t0|1\t1|0",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1\t0|0\t1|0\t0|1",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|1\t0|1\t1|0\t1|1\t0|0\t1|0",
        ];
        let ref_recs: Vec<Rc<dyn RefGTRec>> = ref_lines
            .iter()
            .map(|l| {
                Rc::from(allele_ref_gt_rec_from_parser(&VcfRecGTParser::new(
                    &hr, l, &mp,
                )))
            })
            .collect();
        let ref_gt = RefGT::from_recs(ref_recs);
        let targ_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|0\t1|0\t0|1",
        ];
        let targ_recs: Vec<Rc<dyn GTRec>> = targ_lines
            .iter()
            .map(|l| {
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(&ht, l, &mp))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(targ_recs);
        let indices = MarkerIndices::from_in_targ(&[true, false, false, true], 0, 4);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        let window = Window::new(gen_map.clone(), 1, true, indices, Some(ref_gt), targ);
        let phased: Rc<dyn GT> = Rc::new(window.targ_gt().clone());
        ImpData::new(&par, &window, phased, gen_map.as_ref())
    }

    #[test]
    fn writes_imputed_vcf_records() {
        let imp = imp_data();
        let n_ref = imp.ref_gt().n_markers();
        let n_clusters = imp.n_clusters();
        let state_probs = ImpLS::state_probs(&imp);

        // emit records per target cluster, covering all reference markers like Main does
        let mut out: Vec<u8> = Vec::new();
        for c in 0..n_clusters {
            let mut w = ImputedVcfWriter::new(&imp, 0, n_ref, c);
            w.append_records(&state_probs, &mut out);
        }
        let text = String::from_utf8(out).unwrap();
        let data_lines: Vec<&str> = text.lines().filter(|l| l.starts_with("chr1\t")).collect();
        assert!(!data_lines.is_empty());
        for l in &data_lines {
            assert!(l.contains("GT:DS"), "missing GT:DS in {l}");
            // 3 samples -> 3 GT:DS sample fields after the 9 fixed columns
            assert_eq!(l.split('\t').count(), 9 + 3);
        }
    }
}
