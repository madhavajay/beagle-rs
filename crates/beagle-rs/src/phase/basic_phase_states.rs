//! Port of `phase/BasicPhaseStates.java` — builds the Li & Stephens HMM state alleles for a
//! target sample by assembling up to `maxStates` composite reference haplotypes from the IBS
//! segments found by `PbwtPhaseIbs`, then recording per-cluster allele mismatches.

use crate::beagleutil::CompHapSegment;
use crate::blbutil::{BitArray, Utilities};
use crate::ints::IntIntMap;
use crate::jdk::{PriorityQueue, Random};
use crate::vcf::GT;

use super::{MarkerCluster, PbwtPhaseIbs};

const NIL: i32 = -103;

/// Port of `phase/BasicPhaseStates.java`.
pub struct BasicPhaseStates<'a> {
    ibs_haps: &'a PbwtPhaseIbs<'a>,
    max_states: i32,
    min_steps: i32,
    n_markers: i32,
    hap_to_last_ibs_step: IntIntMap,
    q: PriorityQueue<CompHapSegment>,
    comp_haps: Vec<BitArray>,
}

impl<'a> BasicPhaseStates<'a> {
    /// `new BasicPhaseStates(PbwtPhaseIbs ibsHaps, int maxStates)`.
    pub fn new(ibs_haps: &'a PbwtPhaseIbs<'a>, max_states: i32) -> Self {
        assert!(max_states >= 1, "{max_states}");
        let phase_data = ibs_haps.phase_data();
        let all_haps = ibs_haps.all_haps();
        let markers = all_haps.markers();
        let n_markers = all_haps.n_markers();
        let phase_step = phase_data.fpd().ibs_step();
        // 200 steps and 1 cM
        let min_steps = 200.max((1.0f32 / phase_step).ceil() as i32);
        let n_bits = markers.sum_hap_bits_total();
        let comp_haps: Vec<BitArray> = (0..max_states).map(|_| BitArray::new(n_bits)).collect();
        BasicPhaseStates {
            ibs_haps,
            max_states,
            min_steps,
            n_markers,
            hap_to_last_ibs_step: IntIntMap::new(max_states),
            q: PriorityQueue::new(max_states as usize),
            comp_haps,
        }
    }

    /// `nTargSamples()`.
    pub fn n_targ_samples(&self) -> i32 {
        self.ibs_haps.phase_data().fpd().targ_gt().n_samples()
    }

    /// `nMarkers()`.
    pub fn n_markers(&self) -> i32 {
        self.ibs_haps.phase_data().fpd().targ_gt().n_markers()
    }

    /// `maxStates()`.
    pub fn max_states(&self) -> i32 {
        self.max_states
    }

    /// `ibsStates(MarkerCluster mc, List<int[]> refAtMissingGT, byte[][][] nMismatches)`.
    /// `n_mismatches` has 3 layers `[layer][cluster][state]`.
    pub fn ibs_states(
        &mut self,
        mc: &MarkerCluster,
        ref_at_missing_gt: &mut [Vec<i32>],
        n_mismatches: &mut [Vec<Vec<u8>>],
    ) -> i32 {
        let n_comp_haps = self.set_comp_ref_haps(mc.sample_phase().sample());
        self.copy_data_cluster(mc, n_comp_haps, ref_at_missing_gt, n_mismatches);
        n_comp_haps
    }

    /// `ibsStates(int sample, byte[][][] nMismatches)`. `n_mismatches` has 2 layers
    /// `[layer][marker][state]`.
    pub fn ibs_states_for_sample(&mut self, sample: i32, n_mismatches: &mut [Vec<Vec<u8>>]) -> i32 {
        let n_comp_haps = self.set_comp_ref_haps(sample);
        self.copy_data_sample(sample, n_comp_haps, n_mismatches);
        n_comp_haps
    }

    fn set_comp_ref_haps(&mut self, sample: i32) -> i32 {
        let ibs_haps = self.ibs_haps;
        let steps = ibs_haps.phase_data().fpd().stage1_steps();
        let h1 = sample << 1;
        let h2 = h1 | 0b1;
        self.q.clear();
        self.hap_to_last_ibs_step.clear();
        let n = steps.size();
        for step in 0..n {
            let ibs_hap1 = ibs_haps.ibs_hap(h1, step);
            if ibs_hap1 >= 0 {
                self.add_ibs_hap(ibs_hap1, step);
            }
            let ibs_hap2 = ibs_haps.ibs_hap(h2, step);
            if ibs_hap2 >= 0 {
                self.add_ibs_hap(ibs_hap2, step);
            }
        }
        if self.q.is_empty() {
            self.fill_q_with_random_haps(sample);
        }
        self.copy_final_ref_segs()
    }

    fn add_ibs_hap(&mut self, ibs_hap: i32, step: i32) {
        if self.hap_to_last_ibs_step.get(ibs_hap, NIL) == NIL {
            // hap not currently in q
            self.update_head_of_q();
            let backoff = !self.q.is_empty()
                && step - self.q.peek().expect("nonempty").last_ibs_step() >= self.min_steps;
            if self.q.size() == self.max_states || backoff {
                let mut head = self.q.poll().expect("nonempty");
                let index = head.comp_hap_index();
                let prev_hap = head.hap();
                let prev_start = head.start_marker();
                let next_start = self
                    .ibs_haps
                    .phase_data()
                    .fpd()
                    .stage1_steps()
                    .start((head.last_ibs_step() + step) >> 1);
                self.hap_to_last_ibs_step.remove(head.hap());
                self.ibs_haps.all_haps().copy_to(
                    prev_hap,
                    prev_start,
                    next_start,
                    &mut self.comp_haps[index as usize],
                );
                head.update_segment(ibs_hap, next_start, step);
                self.q.offer(head);
            } else {
                let index = self.q.size();
                let start = 0;
                self.q
                    .offer(CompHapSegment::new(ibs_hap, start, step, index));
            }
        }
        self.hap_to_last_ibs_step.put(ibs_hap, step);
    }

    fn update_head_of_q(&mut self) {
        if self.q.is_empty() {
            return;
        }
        loop {
            let (head_hap, head_last) = {
                let head = self.q.peek().expect("nonempty");
                (head.hap(), head.last_ibs_step())
            };
            let last_ibs_step = self.hap_to_last_ibs_step.get(head_hap, NIL);
            if head_last == last_ibs_step {
                break;
            }
            let mut head = self.q.poll().expect("nonempty");
            head.set_last_ibs_step(last_ibs_step);
            self.q.offer(head);
        }
    }

    fn copy_final_ref_segs(&mut self) -> i32 {
        let ibs_haps = self.ibs_haps;
        let n_markers = self.n_markers;
        let n_comp_haps = self.q.size();
        while let Some(head) = self.q.poll() {
            let index = head.comp_hap_index();
            let hap = head.hap();
            let start_marker = head.start_marker();
            ibs_haps.all_haps().copy_to(
                hap,
                start_marker,
                n_markers,
                &mut self.comp_haps[index as usize],
            );
        }
        n_comp_haps
    }

    // Indices `c`/`j` address parallel [layer][cluster][state] arrays; iterator form obscures
    // the layout and the comp_haps cross-indexing.
    #[allow(clippy::needless_range_loop)]
    fn copy_data_cluster(
        &self,
        mc: &MarkerCluster,
        n_comp_haps: i32,
        ref_at_missing_gt: &mut [Vec<i32>],
        n_mismatches: &mut [Vec<Vec<u8>>],
    ) {
        let markers = self.ibs_haps.all_haps().markers();
        let phase = mc.sample_phase();
        let hap1 = phase.hap1();
        let hap2 = phase.hap2();
        let nc = n_comp_haps as usize;
        let mut miss_index = 0;
        let n_clusters = mc.n_clusters();
        for c in 0..n_clusters {
            let cu = c as usize;
            for k in 0..nc {
                n_mismatches[0][cu][k] = 0;
                n_mismatches[1][cu][k] = 0;
                n_mismatches[2][cu][k] = 0;
            }
            let m_start = mc.cluster_start(c);
            let m_end = mc.cluster_end(c);
            if mc.is_missing_gt_or_masked_het(c) {
                debug_assert_eq!(m_end - m_start, 1);
                let ref_alleles = &mut ref_at_missing_gt[miss_index];
                miss_index += 1;
                for (j, slot) in ref_alleles.iter_mut().enumerate().take(nc) {
                    *slot = markers.allele(&self.comp_haps[j], m_start);
                }
            } else {
                let b_start = markers.sum_hap_bits(m_start);
                let b_end = markers.sum_hap_bits(m_end);
                if hap1.equal(&hap2, b_start, b_end) {
                    for j in 0..nc {
                        if !hap1.equal(&self.comp_haps[j], b_start, b_end) {
                            n_mismatches[0][cu][j] = 1;
                            n_mismatches[1][cu][j] = 1;
                            n_mismatches[2][cu][j] = 1;
                        }
                    }
                } else {
                    // cluster contains a heterozygote genotype
                    for j in 0..nc {
                        if !hap1.equal(&self.comp_haps[j], b_start, b_end) {
                            n_mismatches[1][cu][j] = 1;
                        }
                        if !hap2.equal(&self.comp_haps[j], b_start, b_end) {
                            n_mismatches[2][cu][j] = 1;
                        }
                    }
                }
            }
        }
    }

    // Indices `m`/`j` address parallel [layer][marker][state] arrays plus comp_haps.
    #[allow(clippy::needless_range_loop)]
    fn copy_data_sample(&self, sample: i32, n_comp_haps: i32, n_mismatches: &mut [Vec<Vec<u8>>]) {
        let all_haps = self.ibs_haps.all_haps();
        let markers = all_haps.markers();
        let nc = n_comp_haps as usize;
        let h1 = sample << 1;
        let h2 = h1 | 0b1;
        for m in 0..self.n_markers {
            let mu = m as usize;
            let a1 = all_haps.allele(m, h1);
            let a2 = all_haps.allele(m, h2);
            for j in 0..nc {
                let ref_allele = markers.allele(&self.comp_haps[j], m);
                n_mismatches[0][mu][j] = u8::from(ref_allele != a1);
                n_mismatches[1][mu][j] = u8::from(ref_allele != a2);
            }
        }
    }

    fn fill_q_with_random_haps(&mut self, sample: i32) {
        debug_assert!(self.q.is_empty());
        let ibs_haps = self.ibs_haps;
        let n_haps = ibs_haps.all_haps().n_haps();
        let n_states = (n_haps - 2).min(self.max_states);
        if n_states <= 0 {
            Utilities::exit("ERROR: there is only one sample");
        } else {
            let mut rand = Random::new(ibs_haps.phase_data().seed() + sample as i64);
            let ibs_step = 0;
            let start_marker = 0;
            let mut comp_hap_index = 0;
            for _ in 0..n_states {
                let mut h = rand.next_int_bound(n_haps);
                while (h >> 1) == sample {
                    h = rand.next_int_bound(n_haps);
                }
                if self.hap_to_last_ibs_step.get(h, NIL) == NIL {
                    self.q.add(CompHapSegment::new(
                        h,
                        start_marker,
                        ibs_step,
                        comp_hap_index,
                    ));
                    comp_hap_index += 1;
                    self.hap_to_last_ibs_step.put(h, start_marker);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, HEADER_PREFIX,
    };
    use std::io::Write;
    use std::rc::Rc;

    use super::super::{FixedPhaseData, PhaseData};

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_bps_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..6 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 6],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let gts = [
            "0|0\t0|1\t1|1\t0|1\t0|0\t1|0",
            "0|1\t0|0\t1|0\t1|1\t0|1\t0|0",
            "1|1\t0|1\t0|0\t0|1\t1|0\t0|1",
            "0|0\t1|1\t0|1\t1|0\t0|0\t1|1",
            "0|1\t1|0\t1|1\t0|0\t0|1\t1|0",
            "1|0\t0|1\t0|0\t1|1\t1|0\t0|1",
        ];
        let recs: Vec<Rc<dyn GTRec>> = gts
            .iter()
            .enumerate()
            .map(|(i, g)| {
                let pos = 1_000_000 + i as i32 * 1_000_000;
                let line = format!("chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t{g}");
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(
                    &h, &line, &mp,
                ))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(recs);
        let n = targ.n_markers();
        let indices = MarkerIndices::from_counts(0, n, n);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        let w = Window::new(gen_map, 1, true, indices, None, targ);
        let ped = Pedigree::new(w.targ_gt().samples().clone(), None);
        FixedPhaseData::new(&par, &ped, &w, None)
    }

    #[test]
    fn builds_states_and_records_mismatches() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let phase_ibs = PbwtPhaseIbs::new(&pd, pd.coded_steps(), false);
        let max_states = 4;
        let mut states = BasicPhaseStates::new(&phase_ibs, max_states);
        let n_clusters_max = pd.fpd().stage1_targ_gt().n_markers();

        // buffers sized generously
        let mut ref_at_missing: Vec<Vec<i32>> = (0..n_clusters_max)
            .map(|_| vec![0i32; max_states as usize])
            .collect();
        let mut mismatches: Vec<Vec<Vec<u8>>> = (0..3)
            .map(|_| {
                (0..n_clusters_max)
                    .map(|_| vec![0u8; max_states as usize])
                    .collect()
            })
            .collect();

        let mc = MarkerCluster::new(&pd, 1);
        let n_states = states.ibs_states(&mc, &mut ref_at_missing, &mut mismatches);
        assert!(n_states >= 1 && n_states <= max_states);
        // mismatch flags are 0/1
        for layer in mismatches.iter() {
            for row in layer.iter().take(mc.n_clusters() as usize) {
                for &flag in row.iter().take(n_states as usize) {
                    assert!(flag <= 1);
                }
            }
        }

        // sample variant
        let mut mismatches2: Vec<Vec<Vec<u8>>> = (0..2)
            .map(|_| {
                (0..n_clusters_max)
                    .map(|_| vec![0u8; max_states as usize])
                    .collect()
            })
            .collect();
        let n2 = states.ibs_states_for_sample(2, &mut mismatches2);
        assert!(n2 >= 1 && n2 <= max_states);
    }
}
