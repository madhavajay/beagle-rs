//! Port of `phase/PhaseData.java` — the per-run phasing state: the current phase estimate
//! (`EstPhase`) plus the iteration-dependent parameters (recombination intensity, transition
//! probabilities, allele-mismatch probability, and the likelihood-ratio threshold).
//!
//! Java marks the mutable fields `volatile` for thread safety; the single-threaded Rust port
//! uses `Cell`/`RefCell` so the parameters can be updated through a shared `&PhaseData`.

use std::cell::{Cell, RefCell};
use std::rc::Rc;

use crate::blbutil::FloatArray;
use crate::main_pkg::{li_stephens_p_mismatch, Par};
use crate::vcf::MarkerMap;

use super::{CodedSteps, EstPhase, FixedPhaseData};

/// Port of the `PhaseData.TrProb` inner class.
struct TrProb {
    recomb_intensity: f32,
    p_recomb: FloatArray,
}

impl TrProb {
    fn new(map: &MarkerMap, recomb_intensity: f32) -> Self {
        TrProb {
            recomb_intensity,
            p_recomb: map.p_recomb(recomb_intensity),
        }
    }
}

/// Port of `phase/PhaseData.java`.
pub struct PhaseData {
    est_phase: EstPhase,
    leave_unph_prop: Vec<f32>,
    seed: i64,
    it: Cell<i32>,
    lr_threshold: Cell<f32>,
    tr_prob: RefCell<TrProb>,
    p_mismatch: Cell<f32>,
}

/// `PhaseData.lrThreshold(Par par, int it)`.
fn lr_threshold(par: &Par, it: i32) -> f32 {
    let n_burnin_its = par.burnin();
    let n_its_m1 = par.iterations() - 1;
    if it < n_burnin_its {
        f32::INFINITY
    } else if it == (n_its_m1 + n_burnin_its) {
        1.0
    } else {
        let last_val = 4.0f64;
        let exp = (n_its_m1 - (it - n_burnin_its)) as f64 / n_its_m1 as f64;
        let base = par.initial_lr() as f64 / last_val;
        (last_val * base.powf(exp)) as f32
    }
}

/// `PhaseData.leaveUnphasedProp(FixedPhaseData, EstPhase)`.
fn leave_unphased_prop(fpd: &FixedPhaseData, est_phase: &EstPhase) -> Vec<f32> {
    let n_iterations = fpd.par().iterations();
    (0..fpd.targ_gt().n_samples())
        .map(|s| {
            let cnt = est_phase.get(s).n_unphased();
            (cnt as f64).powf(-1.0 / n_iterations as f64) as f32
        })
        .collect()
}

impl PhaseData {
    /// `new PhaseData(FixedPhaseData fpd, long seed)`.
    pub fn new(fpd: Rc<FixedPhaseData>, seed: i64) -> Self {
        let est_phase = EstPhase::new(fpd.clone(), seed);
        let leave_unph_prop = leave_unphased_prop(&fpd, &est_phase);
        let it = 0;
        let lr_threshold_v = lr_threshold(fpd.par(), it);
        let recomb_intensity = (0.04f32 * fpd.par().ne()) / fpd.n_haps() as f32;
        let tr_prob = TrProb::new(fpd.stage1_map(), recomb_intensity);
        let p_mismatch = li_stephens_p_mismatch(fpd.n_haps());
        PhaseData {
            est_phase,
            leave_unph_prop,
            seed,
            it: Cell::new(it),
            lr_threshold: Cell::new(lr_threshold_v),
            tr_prob: RefCell::new(tr_prob),
            p_mismatch: Cell::new(p_mismatch),
        }
    }

    /// `recombIntensity()`.
    pub fn recomb_intensity(&self) -> f32 {
        self.tr_prob.borrow().recomb_intensity
    }

    /// `updateRecombIntensity(float)`.
    pub fn update_recomb_intensity(&self, recomb_intensity: f32) {
        assert!(
            recomb_intensity > 0.0 && recomb_intensity.is_finite(),
            "{recomb_intensity}"
        );
        let map_recomb = TrProb::new(self.est_phase.fpd().stage1_map(), recomb_intensity);
        *self.tr_prob.borrow_mut() = map_recomb;
    }

    /// `ne()` — the effective population size.
    pub fn ne(&self) -> i64 {
        self.ne_of(self.tr_prob.borrow().recomb_intensity)
    }

    fn ne_of(&self, recomb_intensity: f32) -> i64 {
        ((25.0f32 * recomb_intensity * self.est_phase.fpd().n_haps() as f32) as f64).ceil() as i64
    }

    /// `pRecomb()` — per-marker probability of transitioning to a random HMM state.
    pub fn p_recomb(&self) -> FloatArray {
        self.tr_prob.borrow().p_recomb.clone()
    }

    /// `pMismatch()`.
    pub fn p_mismatch(&self) -> f32 {
        self.p_mismatch.get()
    }

    /// `updatePMismatch(float)`.
    pub fn update_p_mismatch(&self, p_mismatch: f32) {
        assert!(
            (0.0..=1.0).contains(&p_mismatch) && p_mismatch.is_finite(),
            "{p_mismatch}"
        );
        self.p_mismatch.set(p_mismatch);
    }

    /// `incrementIt()`.
    pub fn increment_it(&self) {
        self.it.set(self.it.get() + 1);
        self.lr_threshold
            .set(lr_threshold(self.est_phase.fpd().par(), self.it.get()));
    }

    /// `it()` — the current iteration (initially 0).
    pub fn it(&self) -> i32 {
        self.it.get()
    }

    /// `advanceToFirstPhasingIt()`.
    pub fn advance_to_first_phasing_it(&self) {
        let n_burnin_its = self.est_phase.fpd().par().burnin();
        if self.it.get() < n_burnin_its {
            self.it.set(n_burnin_its);
            self.lr_threshold
                .set(lr_threshold(self.est_phase.fpd().par(), self.it.get()));
        }
    }

    /// `fpd()` — the input data for phasing that is the same in each iteration.
    pub fn fpd(&self) -> &FixedPhaseData {
        self.est_phase.fpd()
    }

    /// `estPhase()` — the estimated phased genotypes.
    pub fn est_phase(&self) -> &EstPhase {
        &self.est_phase
    }

    /// `codedSteps()` — constructed on demand (it can consume substantial memory).
    pub fn coded_steps(&self) -> CodedSteps {
        CodedSteps::new(&self.est_phase)
    }

    /// `leaveUnphasedProp(int sample)`.
    pub fn leave_unphased_prop(&self, sample: i32) -> f32 {
        self.leave_unph_prop[sample as usize]
    }

    /// `lrThreshold()`.
    pub fn lr_threshold(&self) -> f32 {
        self.lr_threshold.get()
    }

    /// `seed()` — an iteration-dependent seed for pseudorandom number generation.
    pub fn seed(&self) -> i64 {
        self.seed + self.it.get() as i64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::main_pkg::Pedigree;
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, GT, HEADER_PREFIX,
    };
    use std::io::Write;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_pd_gt.vcf");
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
    fn parameters_and_iteration_state() {
        let par_burnin;
        let par_iters;
        let pd = {
            let f = fpd();
            par_burnin = f.par().burnin();
            par_iters = f.par().iterations();
            PhaseData::new(Rc::new(f), 99999)
        };

        // Initial iteration 0 < burnin -> threshold is +inf; seed == base seed.
        assert_eq!(pd.it(), 0);
        assert_eq!(pd.seed(), 99999);
        if par_burnin > 0 {
            assert_eq!(pd.lr_threshold(), f32::INFINITY);
        }

        // recombIntensity = 0.04*ne/nHaps; pRecomb has one entry per marker.
        let ri = pd.recomb_intensity();
        assert!(ri > 0.0 && ri.is_finite());
        assert_eq!(pd.p_recomb().size(), pd.fpd().stage1_map().gen_pos().size());
        assert!(pd.ne() >= 1);

        // incrementIt advances iteration and the seed.
        pd.increment_it();
        assert_eq!(pd.it(), 1);
        assert_eq!(pd.seed(), 100000);

        // advanceToFirstPhasingIt jumps to burnin.
        pd.advance_to_first_phasing_it();
        assert_eq!(pd.it(), par_burnin.max(1));

        // updatePMismatch / updateRecombIntensity round-trip.
        pd.update_p_mismatch(0.25);
        assert_eq!(pd.p_mismatch(), 0.25);
        pd.update_recomb_intensity(ri * 2.0);
        assert_eq!(pd.recomb_intensity(), ri * 2.0);

        let _ = par_iters;
    }

    #[test]
    fn coded_steps_built_on_demand() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let cs = pd.coded_steps();
        assert_eq!(cs.targ_samples().size(), 3);
        assert_eq!(pd.est_phase().fpd().stage1_targ_gt().n_samples(), 3);
    }
}
