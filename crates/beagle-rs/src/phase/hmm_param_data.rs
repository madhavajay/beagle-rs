//! Port of `phase/HmmParamData.java` — generates the data (allele-mismatch counts and
//! switch-probability sums) used to estimate the HMM allele-mismatch and recombination-
//! intensity parameters, by running the haploid Li & Stephens forward-backward HMM over the
//! `BasicPhaseStates` for each haplotype of a sample.

use crate::blbutil::FloatArray;

use super::{BasicPhaseStates, HmmUpdater, ParamEstimates, PbwtPhaseIbs};

/// Port of `phase/HmmParamData.java`.
pub struct HmmParamData<'a> {
    n_targ_samples: i32,
    n_markers: i32,
    gen_dist: FloatArray,
    p_recomb: FloatArray,
    al_match: Vec<Vec<Vec<u8>>>, // [2][nMarkers][maxStates]
    states: BasicPhaseStates<'a>,
    fwd: Vec<f32>,
    bwd: Vec<f32>,
    saved_bwd: Vec<Vec<f32>>,
    em_probs: [f32; 2],
    mismatch_cnt: i32,
    sum_mismatch_prob: f64,
    sum_gen_dist: f64,
    sum_switch_prob: f64,
}

impl<'a> HmmParamData<'a> {
    /// `new HmmParamData(PbwtPhaseIbs phaseIbs)`.
    pub fn new(phase_ibs: &'a PbwtPhaseIbs<'a>) -> Self {
        let phase_data = phase_ibs.phase_data();
        let fpd = phase_data.fpd();
        let max_states = fpd.par().phase_states();
        let n_markers = fpd.stage1_targ_gt().n_markers();
        let n_targ_samples = fpd.stage1_targ_gt().n_samples();
        let gen_dist = fpd.stage1_map().gen_dist().clone();
        let p_recomb = phase_data.p_recomb();
        let states = BasicPhaseStates::new(phase_ibs, max_states);
        let p_mismatch = phase_data.p_mismatch();
        let ms = max_states as usize;
        let nm = n_markers as usize;
        HmmParamData {
            n_targ_samples,
            n_markers,
            gen_dist,
            p_recomb,
            al_match: vec![vec![vec![0u8; ms]; nm]; 2],
            states,
            fwd: vec![0.0; ms],
            bwd: vec![0.0; ms],
            saved_bwd: vec![vec![0.0; ms]; nm],
            em_probs: [1.0 - p_mismatch, p_mismatch],
            mismatch_cnt: 0,
            sum_mismatch_prob: 0.0,
            sum_gen_dist: 0.0,
            sum_switch_prob: 0.0,
        }
    }

    /// `nTargSamples()`.
    pub fn n_targ_samples(&self) -> i32 {
        self.n_targ_samples
    }

    /// `addEstimationData(ParamEstimates paramEst)` — flushes the accumulated data and resets.
    pub fn add_estimation_data(&mut self, param_est: &ParamEstimates) {
        param_est.add_mismatch_data(self.mismatch_cnt, self.sum_mismatch_prob);
        param_est.add_switch_data(self.sum_gen_dist, self.sum_switch_prob);
        self.mismatch_cnt = 0;
        self.sum_mismatch_prob = 0.0;
        self.sum_gen_dist = 0.0;
        self.sum_switch_prob = 0.0;
    }

    /// `sumSwitchProbs()`.
    pub fn sum_switch_probs(&self) -> f64 {
        self.sum_switch_prob
    }

    /// `update(int sample)`.
    pub fn update(&mut self, sample: i32) {
        let n_states = self
            .states
            .ibs_states_for_sample(sample, &mut self.al_match);
        if n_states > 1 {
            // otherwise hFactor = nStates/(nStates - 1) is NaN
            self.get_param_data(0, n_states);
            self.get_param_data(1, n_states);
        }
    }

    fn get_param_data(&mut self, layer: usize, n_states: i32) {
        let ns = n_states as usize;
        self.bwd[..ns].fill(1.0);
        let last = (self.n_markers - 1) as usize;
        self.saved_bwd[last][..ns].fill(1.0);
        for m in (0..self.n_markers - 1).rev() {
            let m_p1 = (m + 1) as usize;
            let p_rec = self.p_recomb.get(m + 1);
            HmmUpdater::bwd_update(
                &mut self.bwd,
                p_rec,
                &self.em_probs,
                &self.al_match[layer][m_p1],
                n_states,
            );
            self.saved_bwd[m as usize][..ns].copy_from_slice(&self.bwd[..ns]);
        }
        let h_factor = n_states as f32 / (n_states as f32 - 1.0);
        self.fwd[..ns].fill(1.0 / n_states as f32);
        let mut sum = 1.0f32;
        for m in 0..self.n_markers {
            sum = self.fwd_update(m, layer, n_states, sum, h_factor);
        }
    }

    fn fwd_update(
        &mut self,
        m: i32,
        layer: usize,
        n_states: i32,
        last_sum: f32,
        h_factor: f32,
    ) -> f32 {
        let ns = n_states as usize;
        let mu = m as usize;
        let p_switch = self.p_recomb.get(m);
        let shift = p_switch / n_states as f32;
        let scale = (1.0 - p_switch) / last_sum;
        let no_switch_scale = ((1.0 - p_switch) + shift) / last_sum;
        let mut joint_state_sum = 0.0f32;
        let mut state_sum = 0.0f32;
        let mut fwd_sum = 0.0f32;
        let mut mismatch_sum = 0.0f32;
        for k in 0..ns {
            let al = self.al_match[layer][mu][k];
            let em = self.em_probs[al as usize];
            let bwd_m_k = self.saved_bwd[mu][k];
            joint_state_sum += bwd_m_k * em * no_switch_scale * self.fwd[k];
            self.fwd[k] = em * (scale * self.fwd[k] + shift);
            fwd_sum += self.fwd[k];
            let state_prob = self.fwd[k] * bwd_m_k;
            state_sum += state_prob;
            if al > 0 {
                mismatch_sum += state_prob;
            }
        }
        self.mismatch_cnt += 1;
        self.sum_mismatch_prob += (mismatch_sum / state_sum) as f64;
        let switch_prob = (h_factor * (1.0f32 - joint_state_sum / state_sum)) as f64;
        if switch_prob > 0.0 {
            self.sum_gen_dist += self.gen_dist.get(m) as f64;
            self.sum_switch_prob += switch_prob;
        }
        fwd_sum
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

    use super::super::{FixedPhaseData, PhaseData};

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_hpd_gt.vcf");
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
    fn accumulates_and_feeds_param_estimates() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let phase_ibs = PbwtPhaseIbs::new(&pd, pd.coded_steps(), false);
        let mut hpd = HmmParamData::new(&phase_ibs);
        assert_eq!(hpd.n_targ_samples(), 6);
        for s in 0..hpd.n_targ_samples() {
            hpd.update(s);
        }
        assert!(hpd.sum_switch_probs() >= 0.0);
        let pe = ParamEstimates::new();
        hpd.add_estimation_data(&pe);
        // after flushing, the accumulators reset
        assert_eq!(hpd.sum_switch_probs(), 0.0);
        // estimates are finite/derived without panic
        let _ = pe.p_mismatch();
        let _ = pe.recomb_intensity();
    }
}
