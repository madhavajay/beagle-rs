//! Port of `phase/Stage2Haps.java` — accumulates the stage-2 phased genotypes. High-frequency
//! markers reuse the stage-1 records; low-frequency markers are rebuilt from per-allele rare-
//! carrier haplotype lists populated by `Stage2Baum::set_phased_gt`.

use std::rc::Rc;

use crate::ints::IntArray;
use crate::vcf::{
    allele_ref_gt_rec_from_components, BasicGT, CarrierList, GTRec, Markers, RefGTRec, Samples,
};

use super::FixedPhaseData;

/// Port of `phase/Stage2Haps.java`.
pub struct Stage2Haps {
    fpd: Rc<FixedPhaseData>,
    stage1_recs: Vec<Rc<dyn GTRec>>,
    markers: Markers,
    targ_samples: Samples,
    rare_carriers: Vec<Option<Vec<i32>>>, // [allele offset] -> rare-allele carrier haps, or None
}

fn init_rare_carriers(fpd: &FixedPhaseData) -> Vec<Option<Vec<i32>>> {
    let markers = fpd.targ_gt().markers();
    let size = markers.sum_alleles_total();
    let mut rare: Vec<Option<Vec<i32>>> = vec![None; size as usize];
    let n = fpd.stage1_to2().size();
    for index in 0..=n {
        init_list(fpd, index, &mut rare);
    }
    rare
}

fn init_list(fpd: &FixedPhaseData, index: i32, rare: &mut [Option<Vec<i32>>]) {
    let markers = fpd.targ_gt().markers();
    let stage1_to2 = fpd.stage1_to2();
    let start = if index == 0 {
        0
    } else {
        stage1_to2.get(index - 1) + 1
    };
    let end = if index == stage1_to2.size() {
        markers.size()
    } else {
        stage1_to2.get(index)
    };
    for m in start..end {
        let n_alleles = markers.marker(m).n_alleles();
        let offset = markers.sum_alleles(m);
        for al in 0..n_alleles {
            let ia = fpd.carriers(m, al);
            if !matches!(ia, CarrierList::HighFreq) {
                rare[(offset + al) as usize] = Some(Vec::with_capacity(ia.size() as usize));
            }
        }
    }
}

fn set_major_allele_to_null(hap_indices: &mut [Option<Vec<i32>>]) {
    let mut maj_allele = 0usize;
    for j in 1..hap_indices.len() {
        let len_j = hap_indices[j].as_ref().map_or(0, |v| v.len());
        let len_maj = hap_indices[maj_allele].as_ref().map_or(0, |v| v.len());
        if len_j > len_maj {
            maj_allele = j;
        }
    }
    hap_indices[maj_allele] = None;
}

impl Stage2Haps {
    /// `new Stage2Haps(PhaseData phaseData)`.
    pub fn new(phase_data: &super::PhaseData) -> Self {
        let fpd = phase_data.est_phase().fpd_rc().clone();
        let stage1_recs: Vec<Rc<dyn GTRec>> = phase_data
            .est_phase()
            .to_gt_recs()
            .into_iter()
            .map(|r| Rc::new(r) as Rc<dyn GTRec>)
            .collect();
        let markers = fpd.targ_gt().markers().clone();
        let targ_samples = fpd.targ_gt().samples().clone();
        let rare_carriers = init_rare_carriers(&fpd);
        Stage2Haps {
            fpd,
            stage1_recs,
            markers,
            targ_samples,
            rare_carriers,
        }
    }

    /// `fpd()`.
    pub fn fpd(&self) -> &FixedPhaseData {
        &self.fpd
    }

    /// `setPhasedGT(int marker, int sample, int a1, int a2)`.
    pub fn set_phased_gt(&mut self, marker: i32, sample: i32, a1: i32, a2: i32) {
        let offset = self.markers.sum_alleles(marker);
        if let Some(list1) = self.rare_carriers[(offset + a1) as usize].as_mut() {
            list1.push(sample << 1);
        }
        if let Some(list2) = self.rare_carriers[(offset + a2) as usize].as_mut() {
            list2.push((sample << 1) | 0b1);
        }
    }

    /// `toBasicGT(int start, int end)`.
    pub fn to_basic_gt(&self, start: i32, end: i32) -> BasicGT {
        let recs = self.to_gt_recs(start, end);
        if start == 0 && end == self.markers.size() {
            BasicGT::from_markers_samples(self.markers.clone(), self.targ_samples.clone(), recs)
        } else {
            BasicGT::from_samples(self.targ_samples.clone(), recs)
        }
    }

    /// `toGTRecs(int start, int end)`.
    pub fn to_gt_recs(&self, start: i32, end: i32) -> Vec<Rc<dyn GTRec>> {
        (start..end).map(|m| self.gt_rec(m)).collect()
    }

    fn gt_rec(&self, m: i32) -> Rc<dyn GTRec> {
        let m1 = self.fpd.prev_stage1_marker(m);
        let m2 = self.fpd.stage1_to2().get(m1);
        if m == m2 {
            self.stage1_recs[m1 as usize].clone()
        } else {
            self.stage2_rec(m)
        }
    }

    fn stage2_rec(&self, m: i32) -> Rc<dyn GTRec> {
        let al_start = self.markers.sum_alleles(m);
        let al_end = self.markers.sum_alleles(m + 1);
        let mut hap_indices: Vec<Option<Vec<i32>>> =
            Vec::with_capacity((al_end - al_start) as usize);
        let mut major_allele = -1i32;
        for j in 0..(al_end - al_start) {
            if let Some(list) = &self.rare_carriers[(al_start + j) as usize] {
                let mut arr = list.clone();
                arr.sort_unstable();
                hap_indices.push(Some(arr));
            } else {
                major_allele = j;
                hap_indices.push(None);
            }
        }
        if major_allele == -1 {
            // can occur if all alleles are rare due to a high missing rate
            set_major_allele_to_null(&mut hap_indices);
        }
        let rec = allele_ref_gt_rec_from_components(
            self.markers.marker(m).clone(),
            self.targ_samples.clone(),
            hap_indices,
        );
        let rec: Rc<dyn RefGTRec> = Rc::from(rec);
        rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap, VcfHeader,
        VcfRecGTParser, Window, GT, HEADER_PREFIX,
    };
    use std::io::Write;

    use super::super::PhaseData;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_s2h_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..4 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 4],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let recs: Vec<Rc<dyn GTRec>> = (0..4)
            .map(|i| {
                let pos = 1_000_000 + i * 1_000_000;
                let line = format!("chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|0");
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
    fn full_range_basic_gt_has_expected_shape() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let n_markers = pd.fpd().targ_gt().n_markers();
        let s2h = Stage2Haps::new(&pd);
        assert_eq!(s2h.fpd().targ_gt().n_samples(), 4);
        let gt = s2h.to_basic_gt(0, n_markers);
        assert_eq!(gt.n_markers(), n_markers);
        assert_eq!(gt.n_samples(), 4);
        // every stage-1 marker maps to a stage1 record (here all markers are stage-1 markers)
        let recs = s2h.to_gt_recs(0, n_markers);
        assert_eq!(recs.len() as i32, n_markers);
    }
}
