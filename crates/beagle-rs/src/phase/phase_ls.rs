//! Port of `phase/PhaseLS.java` — the static driver that updates genotype phase with the
//! haploid Li & Stephens HMM. `runStage1` re-phases the high-frequency markers (optionally
//! EM-estimating the HMM parameters during burn-in); `runStage2` phases all markers.
//!
//! Java fans the per-sample work across `nthreads` executor threads sharing an atomic counter;
//! the work is deterministic per sample (seeded from `seed + sample`), so the single-threaded
//! port loops over the samples in order with one HMM instance.

use crate::blbutil::Utilities;
use crate::jdk::Random;

use super::{
    HmmParamData, LowFreqPhaseIbs, ParamEstimates, PbwtPhaseIbs, PhaseBaum, PhaseBaum2, PhaseData,
    Stage2Baum, Stage2Haps,
};

/// Port of `phase/PhaseLS.java`.
pub struct PhaseLS;

impl PhaseLS {
    /// `runStage1(PhaseData pd)` — updates genotype phase estimates at the stage-1 markers.
    pub fn run_stage1(pd: &PhaseData) {
        let it = pd.it();
        let n_samples = pd.fpd().targ_gt().n_samples();
        let n_burnin_its = pd.fpd().par().burnin();
        let phase_ibs = pbwt_phase_ibs(pd);
        if pd.fpd().par().em() {
            let mut rand = Random::new(pd.seed());
            if it == 0 {
                initialize_parameters(&phase_ibs, &mut rand);
            } else if it < n_burnin_its {
                update_parameters(&phase_ibs, &mut rand);
            }
        }
        let mut baum = PhaseBaum2::new(&phase_ibs);
        for s in 0..n_samples {
            baum.phase(s);
        }
    }

    /// `runStage2(PhaseData pd)` — returns phased genotypes at all markers.
    pub fn run_stage2(pd: &PhaseData) -> Stage2Haps {
        let n_samples = pd.fpd().targ_gt().n_samples();
        let phase_ibs = LowFreqPhaseIbs::new(pd);
        let mut stage2_haps = Stage2Haps::new(pd);
        {
            let mut baum = Stage2Baum::new(&phase_ibs, &mut stage2_haps);
            for s in 0..n_samples {
                baum.phase(s);
            }
        }
        stage2_haps
    }
}

fn pbwt_phase_ibs(pd: &PhaseData) -> PbwtPhaseIbs<'_> {
    let use_bwd = (pd.it() & 1) == 0;
    PbwtPhaseIbs::new(pd, pd.coded_steps(), use_bwd)
}

fn initialize_parameters(phase_ibs: &PbwtPhaseIbs, rand: &mut Random) {
    let pd = phase_ibs.phase_data();
    let mut prev_rec_int = pd.recomb_intensity();
    let max_initial_its = 15;
    for _ in 0..max_initial_its {
        update_parameters(phase_ibs, rand);
        let rec_int = pd.recomb_intensity();
        // Java: Math.abs(recInt - prevRecInt) <= 0.1*prevRecInt (float compared as double)
        if (rec_int - prev_rec_int).abs() as f64 <= 0.1f64 * prev_rec_int as f64 {
            break;
        }
        prev_rec_int = rec_int;
    }
}

fn update_parameters(phase_ibs: &PbwtPhaseIbs, rand: &mut Random) {
    let param_est = get_param_est(phase_ibs, rand);
    let pd = phase_ibs.phase_data();
    let prev_p_mismatch = pd.p_mismatch();
    let p_mismatch = param_est.p_mismatch();
    let recomb_intensity = param_est.recomb_intensity();
    if p_mismatch.is_finite() && p_mismatch > prev_p_mismatch {
        pd.update_p_mismatch(p_mismatch);
    }
    if recomb_intensity.is_finite() && recomb_intensity > 0.0 {
        pd.update_recomb_intensity(recomb_intensity);
    }
}

fn get_param_est(phase_ibs: &PbwtPhaseIbs, rand: &mut Random) -> ParamEstimates {
    let param_est = ParamEstimates::new();
    let pd = phase_ibs.phase_data();
    let fpd = pd.fpd();
    let sample_indices = samples_to_analyze(pd, rand);
    let n_threads = fpd.par().nthreads().min(sample_indices.len() as i32);
    let max_sum = 20000.0 / n_threads as f64;
    let min_indices = 50;
    let mut hpd = HmmParamData::new(phase_ibs);
    let mut index = 0i32;
    while (hpd.sum_switch_probs() < max_sum || index < min_indices)
        && index < sample_indices.len() as i32
    {
        hpd.update(sample_indices[index as usize]);
        index += 1;
        hpd.add_estimation_data(&param_est);
    }
    param_est
}

fn samples_to_analyze(pd: &PhaseData, rand: &mut Random) -> Vec<i32> {
    let max_samples_to_analyze = 500;
    let n_targ_samples = pd.fpd().targ_gt().n_samples();
    let mut ia: Vec<i32> = (0..n_targ_samples).collect();
    if n_targ_samples <= max_samples_to_analyze {
        ia
    } else {
        Utilities::shuffle_n(&mut ia, max_samples_to_analyze, rand);
        ia.truncate(max_samples_to_analyze as usize);
        ia
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
        let gt = std::env::temp_dir().join("beagle_rs_pls_gt.vcf");
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
            "1|1\t0|1\t.|.\t0|1\t1|0\t0|1",
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
    fn stage1_then_stage2_fully_phases() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        // a couple of phasing iterations
        PhaseLS::run_stage1(&pd);
        pd.increment_it();
        PhaseLS::run_stage1(&pd);

        let n_markers = pd.fpd().targ_gt().n_markers();
        let stage2 = PhaseLS::run_stage2(&pd);
        let gt = stage2.to_basic_gt(0, n_markers);
        assert_eq!(gt.n_markers(), n_markers);
        assert_eq!(gt.n_samples(), 6);
        // all genotypes phased and non-missing (the missing one at sample 2 / marker 2 imputed)
        for m in 0..n_markers {
            for hap in 0..12 {
                assert!(gt.allele(m, hap) >= 0, "m={m} hap={hap} unphased");
            }
        }
    }
}
