//! Port of `phase/LowFreqPhaseIbs.java` — pairs a forward and a backward
//! `LowFreqPbwtPhaseIbs` so callers can query the IBS haplotype from either PBWT direction.

use crate::vcf::XRefGT;

use super::{LowFreqPbwtPhaseIbs, PhaseData};

/// Port of `phase/LowFreqPhaseIbs.java`.
pub struct LowFreqPhaseIbs<'a> {
    phase_data: &'a PhaseData,
    all_haps: XRefGT,
    fwd_phase_ibs: LowFreqPbwtPhaseIbs<'a>,
    bwd_phase_ibs: LowFreqPbwtPhaseIbs<'a>,
}

impl<'a> LowFreqPhaseIbs<'a> {
    /// `new LowFreqPhaseIbs(PhaseData phaseData)`.
    pub fn new(phase_data: &'a PhaseData) -> Self {
        // Java builds the coded steps once and shares them; the build is deterministic, so
        // rebuilding per phaser (each consumes its `CodedSteps`) yields the identical result.
        let coded_steps = phase_data.coded_steps();
        let all_haps = coded_steps.all_haps().clone();
        let fwd_phase_ibs = LowFreqPbwtPhaseIbs::new(phase_data, coded_steps, false);
        let bwd_phase_ibs = LowFreqPbwtPhaseIbs::new(phase_data, phase_data.coded_steps(), true);
        LowFreqPhaseIbs {
            phase_data,
            all_haps,
            fwd_phase_ibs,
            bwd_phase_ibs,
        }
    }

    /// `phaseData()`.
    pub fn phase_data(&self) -> &PhaseData {
        self.phase_data
    }

    /// `allHaps()`.
    pub fn all_haps(&self) -> &XRefGT {
        &self.all_haps
    }

    /// `fwdIbsHap(int targHap, int step)` — IBS hap from the increasing-index PBWT, or `-1`.
    pub fn fwd_ibs_hap(&self, targ_hap: i32, step: i32) -> i32 {
        self.fwd_phase_ibs.ibs_hap(targ_hap, step)
    }

    /// `bwdIbsHap(int targHap, int step)` — IBS hap from the decreasing-index PBWT, or `-1`.
    pub fn bwd_ibs_hap(&self, targ_hap: i32, step: i32) -> i32 {
        self.bwd_phase_ibs.ibs_hap(targ_hap, step)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, GT, HEADER_PREFIX,
    };
    use std::io::Write;
    use std::rc::Rc;

    use super::super::FixedPhaseData;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_lfpi_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..5 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 5],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let gts = [
            "0|0\t0|1\t1|1\t0|1\t0|0",
            "0|1\t0|0\t1|0\t1|1\t0|1",
            "1|1\t0|1\t0|0\t0|1\t1|0",
            "0|0\t1|1\t0|1\t1|0\t0|0",
            "0|1\t1|0\t1|1\t0|0\t0|1",
            "1|0\t0|1\t0|0\t1|1\t1|0",
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
    fn fwd_and_bwd_queries_are_valid() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let ibs = LowFreqPhaseIbs::new(&pd);
        let n_haps = ibs.phase_data().fpd().targ_gt().n_haps();
        let n_steps = pd.fpd().stage1_steps().size();
        assert_eq!(ibs.all_haps().n_haps(), n_haps);
        for step in 0..n_steps {
            for hap in 0..n_haps {
                let f = ibs.fwd_ibs_hap(hap, step);
                let b = ibs.bwd_ibs_hap(hap, step);
                assert!(f == -1 || (0..n_haps).contains(&f));
                assert!(b == -1 || (0..n_haps).contains(&b));
            }
        }
    }
}
