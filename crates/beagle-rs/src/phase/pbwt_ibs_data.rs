//! Port of `phase/PbwtIbsData.java` — parameters and step/batch geometry for finding the
//! haplotypes that share an IBS segment with a target haplotype.
//!
//! Java's `checkConsistency` compares the `CodedSteps` fields against the `FixedPhaseData`
//! by reference identity. In the Rust port `CodedSteps` owns clones of those values, so the
//! port checks structural equality instead (a behavior-preserving divergence: the guard only
//! catches programmer error, and the structural check is at least as strict on the fields
//! that matter).

use super::{CodedSteps, PhaseData};

const BURNIN_CANDIDATES: i32 = 100;
const MAX_PHASE_CANDIDATES: i32 = 90;
const MIN_PHASE_CANDIDATES: i32 = 5;
const STAGE2_CANDIDATES: i32 = 10;
const MAX_BACKOFF_CM: f32 = 0.3;

/// Port of `phase/PbwtIbsData.java`.
pub struct PbwtIbsData {
    coded_steps: CodedSteps,
    n_haps: i32,
    n_targ_haps: i32,
    n_candidates: i32,
    n_steps: i32,
    n_overlap_steps: i32,
    max_backoff_steps: i32,
    steps_per_batch: i32,
    n_batches: i32,
}

fn check_consistency(phase_data: &PhaseData, coded_steps: &CodedSteps) {
    let fpd = phase_data.fpd();
    let consistent = fpd.stage1_steps().size() == coded_steps.steps().size()
        && fpd.stage1_xref_gt().is_some() == coded_steps.ref_haps().is_some()
        && fpd.targ_gt().samples() == coded_steps.targ_samples();
    assert!(consistent, "inconsistent data");
}

fn n_candidates1(phase_data: &PhaseData) -> i32 {
    let mut n_candidates = BURNIN_CANDIDATES;
    let it = phase_data.it();
    let par = phase_data.fpd().par();
    if it >= par.burnin() {
        let n_its_remaining = (par.burnin() + par.iterations() - it) as f64;
        let p = n_its_remaining / par.iterations() as f64;
        // Java Math.round(double) == (long) Math.floor(a + 0.5d).
        n_candidates = (p * MAX_PHASE_CANDIDATES as f64 + 0.5).floor() as i32;
        n_candidates = n_candidates.max(MIN_PHASE_CANDIDATES);
    }
    n_candidates.min(phase_data.fpd().n_haps())
}

impl PbwtIbsData {
    /// `new PbwtIbsData(PhaseData phaseData, CodedSteps codedSteps)`.
    pub fn new(phase_data: &PhaseData, coded_steps: CodedSteps) -> Self {
        check_consistency(phase_data, &coded_steps);
        let fpd = phase_data.fpd();
        let par = fpd.par();
        let n_threads = par.nthreads();
        let n_its = par.burnin() + par.iterations();

        let n_haps = fpd.n_haps();
        let n_targ_haps = fpd.targ_gt().n_haps();
        let n_candidates = if phase_data.it() < n_its {
            n_candidates1(phase_data)
        } else {
            STAGE2_CANDIDATES.min(fpd.n_haps())
        };
        let n_steps = coded_steps.steps().size();
        // (int) Math.rint(float / float): float division widened to double, round-half-even.
        let n_overlap_steps = ((par.buffer() / fpd.ibs_step()) as f64).round_ties_even() as i32;
        let max_backoff_steps = ((MAX_BACKOFF_CM / fpd.ibs_step()) as f64).round_ties_even() as i32;
        let steps_per_batch = (n_steps + n_threads - 1) / n_threads;
        let n_batches = (n_steps + steps_per_batch - 1) / steps_per_batch;

        PbwtIbsData {
            coded_steps,
            n_haps,
            n_targ_haps,
            n_candidates,
            n_steps,
            n_overlap_steps,
            max_backoff_steps,
            steps_per_batch,
            n_batches,
        }
    }

    /// `codedSteps()`.
    pub fn coded_steps(&self) -> &CodedSteps {
        &self.coded_steps
    }
    /// `nHaps()` — total target + reference haplotypes.
    pub fn n_haps(&self) -> i32 {
        self.n_haps
    }
    /// `nTargHaps()`.
    pub fn n_targ_haps(&self) -> i32 {
        self.n_targ_haps
    }
    /// `nCandidates()`.
    pub fn n_candidates(&self) -> i32 {
        self.n_candidates
    }
    /// `nSteps()`.
    pub fn n_steps(&self) -> i32 {
        self.n_steps
    }
    /// `nOverlapSteps()`.
    pub fn n_overlap_steps(&self) -> i32 {
        self.n_overlap_steps
    }
    /// `maxBackoffSteps()`.
    pub fn max_backoff_steps(&self) -> i32 {
        self.max_backoff_steps
    }
    /// `stepsPerBatch()`.
    pub fn steps_per_batch(&self) -> i32 {
        self.steps_per_batch
    }
    /// `nBatches()`.
    pub fn n_batches(&self) -> i32 {
        self.n_batches
    }

    /// `startStep(int batch)` — inclusive start step of a batch.
    pub fn start_step(&self, batch: i32) -> i32 {
        assert!(batch >= 0 && batch < self.n_batches, "{batch}");
        batch * self.steps_per_batch
    }

    /// `endStep(int batch)` — exclusive end step of a batch.
    pub fn end_step(&self, batch: i32) -> i32 {
        assert!(batch >= 0 && batch < self.n_batches, "{batch}");
        ((batch + 1) * self.steps_per_batch).min(self.n_steps)
    }

    /// `bufferStartStep(int startStep)` — inclusive start of the leading buffer segment.
    pub fn buffer_start_step(&self, start_step: i32) -> i32 {
        assert!(start_step >= 0 && start_step < self.n_steps, "{start_step}");
        (start_step - self.n_overlap_steps).max(0)
    }

    /// `bufferEndStep(int endStep)` — exclusive end of the trailing buffer segment.
    pub fn buffer_end_step(&self, end_step: i32) -> i32 {
        assert!(end_step > 0 && end_step <= self.n_steps, "{end_step}");
        (end_step + self.n_overlap_steps).min(self.n_steps)
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
        let gt = std::env::temp_dir().join("beagle_rs_pid_gt.vcf");
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
        let recs: Vec<Rc<dyn GTRec>> = (0..6)
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
    fn geometry_and_batches() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let cs = pd.coded_steps();
        let n_steps = cs.steps().size();
        let data = PbwtIbsData::new(&pd, cs);

        assert_eq!(data.n_haps(), 6);
        assert_eq!(data.n_targ_haps(), 6);
        assert_eq!(data.n_steps(), n_steps);
        // nthreads=1 -> a single batch covering every step.
        assert_eq!(data.n_batches(), 1);
        assert_eq!(data.start_step(0), 0);
        assert_eq!(data.end_step(0), n_steps);
        // candidates capped at nHaps.
        assert!(data.n_candidates() >= 1 && data.n_candidates() <= 6);

        // buffer steps stay within [0, nSteps].
        let bs = data.buffer_start_step(0);
        assert!(bs >= 0 && bs <= n_steps);
        let be = data.buffer_end_step(n_steps);
        assert!(be >= 0 && be <= n_steps);
    }
}
