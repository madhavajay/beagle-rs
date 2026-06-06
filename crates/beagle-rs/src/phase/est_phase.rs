//! Port of `phase/EstPhase.java` — stores the fixed phasing input together with the current
//! estimated phased genotypes for each target sample, updated in place across iterations.
//!
//! Java uses an `AtomicReferenceArray<SamplePhase>` for thread safety. The single-threaded
//! Rust port stores the estimates in a `RefCell<Vec<SamplePhase>>`; `get` returns an owned
//! clone (Java hands out a shared reference to the same object — callers mutate it and call
//! `set` to store it back, which the clone/`set` pair reproduces).

use std::cell::RefCell;
use std::rc::Rc;

use crate::vcf::{BitArrayRefGTRec, XRefGT};

use super::{FixedPhaseData, PbwtPhaser, SamplePhase};

/// Port of `phase/EstPhase.java`.
pub struct EstPhase {
    fpd: Rc<FixedPhaseData>,
    phase: RefCell<Vec<SamplePhase>>,
}

impl EstPhase {
    /// `new EstPhase(FixedPhaseData fpd, long seed)`.
    pub fn new(fpd: Rc<FixedPhaseData>, seed: i64) -> Self {
        let phase = PbwtPhaser::init_phase(&fpd, seed);
        EstPhase {
            fpd,
            phase: RefCell::new(phase),
        }
    }

    /// `fpd()` — the input data for phasing that is the same in each iteration.
    pub fn fpd(&self) -> &FixedPhaseData {
        &self.fpd
    }

    /// Shared handle to the fixed phasing input (Java passes the same `FixedPhaseData`
    /// object to every consumer).
    pub fn fpd_rc(&self) -> &Rc<FixedPhaseData> {
        &self.fpd
    }

    /// `set(int sample, SamplePhase)`.
    pub fn set(&self, sample: i32, sample_phase: SamplePhase) {
        self.phase.borrow_mut()[sample as usize] = sample_phase;
    }

    /// `get(int sample)` — the estimated phase for `sample` (an owned clone).
    pub fn get(&self, sample: i32) -> SamplePhase {
        self.phase.borrow()[sample as usize].clone()
    }

    /// `phasedHaps()` — the current estimate as haplotype-major `XRefGT`.
    pub fn phased_haps(&self) -> XRefGT {
        let phase = self.phase.borrow();
        XRefGT::from(self.fpd.stage1_targ_gt().samples(), &phase)
    }

    /// `toGTRecs()` — the current estimate as row-major per-marker `BitArrayRefGTRec`s.
    pub fn to_gt_recs(&self) -> Vec<BitArrayRefGTRec> {
        BitArrayRefGTRec::to_bit_array_ref_gt_recs(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::IntArray;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, GT, HEADER_PREFIX,
    };
    use std::io::Write;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_est_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..3 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 3],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let recs: Vec<Rc<dyn GTRec>> = (0..4)
            .map(|i| {
                let pos = 1_000_000 + i * 1_000_000;
                let line = format!("chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1");
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
    fn phased_haps_and_gt_recs_round_trip() {
        let est = EstPhase::new(Rc::new(fpd()), 99999);
        let n = est.fpd().stage1_targ_gt().n_markers();

        // phasedHaps(): one XRefGT over the 3 samples, genotypes preserved.
        let haps = est.phased_haps();
        assert_eq!(haps.n_samples(), 3);
        assert_eq!(haps.n_markers(), n);

        // toGTRecs(): one record per marker, row-major; must agree with phasedHaps().
        let recs = est.to_gt_recs();
        assert_eq!(recs.len() as i32, n);
        for m in 0..n {
            let rec = &recs[m as usize];
            assert_eq!(rec.size(), 6); // 3 samples * 2 haps
            for hap in 0..6 {
                assert_eq!(rec.get(hap), haps.allele(m, hap), "m={m} hap={hap}");
            }
            // S0 = 0|0, S2 = 1|1; S1 het in some phase.
            assert_eq!((rec.get(0), rec.get(1)), (0, 0));
            assert_eq!((rec.get(4), rec.get(5)), (1, 1));
            let s1 = (rec.get(2), rec.get(3));
            assert!(s1 == (0, 1) || s1 == (1, 0));
        }
    }

    #[test]
    fn get_set_round_trips_sample_phase() {
        let est = EstPhase::new(Rc::new(fpd()), 99999);
        let sp = est.get(1);
        // get returns an independent clone; mutating it and setting it back updates the store.
        let mut sp2 = est.get(1);
        sp2.set_allele1(0, sp.allele2(0));
        sp2.set_allele2(0, sp.allele1(0));
        est.set(1, sp2);
        let after = est.get(1);
        assert_eq!(after.allele1(0), sp.allele2(0));
        assert_eq!(after.allele2(0), sp.allele1(0));
    }
}
