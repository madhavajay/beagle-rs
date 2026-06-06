//! Port of `phase/HmmStateProbs.java` — runs the haploid Li & Stephens forward-backward HMM
//! over the `LowFreqPhaseStates` for one target haplotype and returns, per marker, the
//! reference haplotype and posterior probability of each hidden state.

use crate::blbutil::FloatArray;

use super::{LowFreqPhaseIbs, LowFreqPhaseStates};

/// Port of `phase/HmmStateProbs.java`.
pub struct HmmStateProbs<'a> {
    states: LowFreqPhaseStates<'a>,
    p_recomb: FloatArray,
    mismatch: Vec<Vec<u8>>,
    bwd: Vec<f32>,
    p_mismatch: [f32; 2],
}

impl<'a> HmmStateProbs<'a> {
    /// `new HmmStateProbs(LowFreqPhaseIbs phaseIbs)`.
    pub fn new(phase_ibs: &'a LowFreqPhaseIbs<'a>) -> Self {
        let phase_data = phase_ibs.phase_data();
        let n_markers = phase_data.fpd().stage1_targ_gt().n_markers();
        let max_states = phase_data.fpd().par().phase_states() / 2;
        let states = LowFreqPhaseStates::new(phase_ibs, max_states);
        let p_recomb = phase_data.p_recomb();
        let p_miss = phase_data.p_mismatch();
        HmmStateProbs {
            states,
            p_recomb,
            mismatch: vec![vec![0u8; max_states as usize]; n_markers as usize],
            bwd: vec![0.0f32; max_states as usize],
            p_mismatch: [1.0 - p_miss, p_miss],
        }
    }

    /// `nMarkers()`.
    pub fn n_markers(&self) -> i32 {
        self.mismatch.len() as i32
    }

    /// `nTargHaps()`.
    pub fn n_targ_haps(&self) -> i32 {
        self.states.n_targ_haps()
    }

    /// `maxStates()`.
    pub fn max_states(&self) -> i32 {
        self.bwd.len() as i32
    }

    /// `run(int targHap, int[][] refHaps, float[][] stateProbs)` — fills `ref_haps[m][j]` and
    /// the posterior `state_probs[m][j]`, returning the number of states per marker.
    pub fn run(
        &mut self,
        targ_hap: i32,
        ref_haps: &mut [Vec<i32>],
        state_probs: &mut [Vec<f32>],
    ) -> i32 {
        let n_states = self
            .states
            .ibs_states(targ_hap, ref_haps, &mut self.mismatch);
        self.run_fwd(state_probs, n_states);
        self.run_bwd(state_probs, n_states);
        n_states
    }

    // `j` indexes the per-state probability/mismatch columns in lockstep.
    #[allow(clippy::needless_range_loop)]
    fn run_fwd(&self, probs: &mut [Vec<f32>], n_states: i32) {
        let ns = n_states as usize;
        let mut last_sum = 0.0f32;
        for j in 0..ns {
            probs[0][j] = self.p_mismatch[self.mismatch[0][j] as usize];
            last_sum += probs[0][j];
        }
        for m in 1..probs.len() {
            let p_rec = self.p_recomb.get(m as i32);
            let shift = p_rec / n_states as f32;
            let scale = (1.0 - p_rec) / last_sum;
            last_sum = 0.0;
            let (head, tail) = probs.split_at_mut(m);
            let prev = &head[m - 1];
            let cur = &mut tail[0];
            for j in 0..ns {
                let em = self.p_mismatch[self.mismatch[m][j] as usize];
                cur[j] = em * (scale * prev[j] + shift);
                last_sum += cur[j];
            }
        }
    }

    #[allow(clippy::needless_range_loop)]
    fn run_bwd(&mut self, probs: &mut [Vec<f32>], n_states: i32) {
        let ns = n_states as usize;
        let incl_end = probs.len() - 1;
        self.bwd[..ns].fill(1.0 / n_states as f32);
        for m in (0..incl_end).rev() {
            let m_p1 = m + 1;
            let mut sum = 0.0f32;
            for j in 0..ns {
                self.bwd[j] *= self.p_mismatch[self.mismatch[m_p1][j] as usize];
                sum += self.bwd[j];
            }
            let p_rec = self.p_recomb.get(m_p1 as i32);
            let scale = (1.0 - p_rec) / sum;
            let shift = p_rec / n_states as f32;
            sum = 0.0;
            for j in 0..ns {
                self.bwd[j] = scale * self.bwd[j] + shift;
                probs[m][j] *= self.bwd[j];
                sum += probs[m][j];
            }
            for j in 0..ns {
                probs[m][j] /= sum;
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
        VcfHeader, VcfRecGTParser, Window, GT, HEADER_PREFIX,
    };
    use std::io::Write;
    use std::rc::Rc;

    use super::super::{FixedPhaseData, PhaseData};

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_hsp_gt.vcf");
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
    fn posterior_probs_sum_to_one_per_marker() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let ibs = LowFreqPhaseIbs::new(&pd);
        let mut hsp = HmmStateProbs::new(&ibs);
        let n_markers = hsp.n_markers();
        let max_states = hsp.max_states();
        assert_eq!(n_markers, pd.fpd().stage1_targ_gt().n_markers());

        let mut ref_haps: Vec<Vec<i32>> = (0..n_markers)
            .map(|_| vec![0i32; max_states as usize])
            .collect();
        let mut probs: Vec<Vec<f32>> = (0..n_markers)
            .map(|_| vec![0.0f32; max_states as usize])
            .collect();

        let n_states = hsp.run(0, &mut ref_haps, &mut probs);
        assert!(n_states >= 1 && n_states <= max_states);
        // backward pass normalizes every marker except the last to sum 1
        for (m, row) in probs.iter().enumerate().take((n_markers - 1) as usize) {
            let s: f32 = row.iter().take(n_states as usize).sum();
            assert!((s - 1.0).abs() < 1e-3, "marker {m} sum {s}");
        }
    }
}
