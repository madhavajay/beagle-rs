//! Port of `imp/StateProbsFactory.java` — builds a [`StateProbs`] for a target haplotype,
//! keeping only the HMM states whose probability at a marker or the following marker exceeds a
//! threshold (so missing-marker probabilities can later be linearly interpolated).

use super::StateProbs;

/// Port of `imp/StateProbsFactory.java`.
pub struct StateProbsFactory {
    n_targ_markers: i32,
    n_targ_markers_m1: i32,
}

impl StateProbsFactory {
    /// `new StateProbsFactory(int nTargMarkers)`.
    pub fn new(n_targ_markers: i32) -> Self {
        assert!(n_targ_markers > 0, "{n_targ_markers}");
        StateProbsFactory {
            n_targ_markers,
            n_targ_markers_m1: n_targ_markers - 1,
        }
    }

    /// `nTargMarkers()`.
    pub fn n_targ_markers(&self) -> i32 {
        self.n_targ_markers
    }

    fn threshold(n_states: i32) -> f32 {
        0.005f32.min(0.9999f32 / n_states as f32)
    }

    /// `stateProbs(int targHap, int nStates, int[][] hapIndices, float[][] stateProbs)`.
    pub fn state_probs(
        &self,
        targ_hap: i32,
        n_states: i32,
        hap_indices: &[Vec<i32>],
        state_probs: &[Vec<f32>],
    ) -> BasicStateProbs {
        assert!(
            hap_indices.len() as i32 == self.n_targ_markers,
            "{}",
            hap_indices.len()
        );
        assert!(
            state_probs.len() as i32 == self.n_targ_markers,
            "{}",
            state_probs.len()
        );
        let threshold = Self::threshold(n_states);
        let nm = self.n_targ_markers as usize;
        let mut haps: Vec<Vec<i32>> = Vec::with_capacity(nm);
        let mut probs: Vec<Vec<f32>> = Vec::with_capacity(nm);
        let mut probs_p1: Vec<Vec<f32>> = Vec::with_capacity(nm);
        for m in 0..nm {
            let m_p1 = if (m as i32) < self.n_targ_markers_m1 {
                m + 1
            } else {
                m
            };
            let mut hl = Vec::new();
            let mut pl = Vec::new();
            let mut p1l = Vec::new();
            for j in 0..n_states as usize {
                if state_probs[m][j] > threshold || state_probs[m_p1][j] > threshold {
                    hl.push(hap_indices[m][j]);
                    pl.push(state_probs[m][j]);
                    p1l.push(state_probs[m_p1][j]);
                }
            }
            haps.push(hl);
            probs.push(pl);
            probs_p1.push(p1l);
        }
        BasicStateProbs {
            targ_hap,
            haps,
            probs,
            probs_p1,
        }
    }
}

/// Port of the `StateProbsFactory.BasicStateProbs` inner class.
pub struct BasicStateProbs {
    targ_hap: i32,
    haps: Vec<Vec<i32>>,
    probs: Vec<Vec<f32>>,
    probs_p1: Vec<Vec<f32>>,
}

impl StateProbs for BasicStateProbs {
    fn targ_hap(&self) -> i32 {
        self.targ_hap
    }

    fn n_targ_markers(&self) -> i32 {
        self.haps.len() as i32
    }

    fn n_states(&self, marker: i32) -> i32 {
        self.haps[marker as usize].len() as i32
    }

    fn ref_hap(&self, marker: i32, index: i32) -> i32 {
        self.haps[marker as usize][index as usize]
    }

    fn probs(&self, marker: i32, index: i32) -> f32 {
        self.probs[marker as usize][index as usize]
    }

    fn probs_p1(&self, marker: i32, index: i32) -> f32 {
        self.probs_p1[marker as usize][index as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_states_above_threshold() {
        let n_markers = 3;
        let n_states = 4;
        let factory = StateProbsFactory::new(n_markers);
        // threshold = min(0.005, 0.9999/4) = 0.005
        let hap_indices: Vec<Vec<i32>> = vec![
            vec![10, 11, 12, 13],
            vec![20, 21, 22, 23],
            vec![30, 31, 32, 33],
        ];
        let state_probs: Vec<Vec<f32>> = vec![
            vec![0.9, 0.001, 0.001, 0.098],
            vec![0.5, 0.5, 0.0, 0.0],
            vec![0.001, 0.001, 0.001, 0.997],
        ];
        let sp = factory.state_probs(0, n_states, &hap_indices, &state_probs);
        assert_eq!(sp.targ_hap(), 0);
        assert_eq!(sp.n_targ_markers(), 3);
        // marker 0: states above threshold at m or m+1. state1/state2 are 0.001 at m0 and
        // 0.5/0.0 at m1 -> state1 kept (0.5 at m1), state2 dropped (0.001 & 0.0).
        let m0_haps: Vec<i32> = (0..sp.n_states(0)).map(|i| sp.ref_hap(0, i)).collect();
        assert!(m0_haps.contains(&10)); // 0.9
        assert!(m0_haps.contains(&13)); // 0.098
        assert!(m0_haps.contains(&11)); // 0.5 at m+1
        assert!(!m0_haps.contains(&12)); // below threshold at m and m+1
                                         // probsP1 at the last marker equals probs (mP1 == m)
        for i in 0..sp.n_states(2) {
            assert_eq!(sp.probs(2, i), sp.probs_p1(2, i));
        }
    }
}
