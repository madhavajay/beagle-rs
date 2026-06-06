//! Port of `phase/HmmUpdater.java` — single-marker forward/backward HMM updates.
//!
//! All arithmetic is `f32` in the same operation order as Java so results are bit-identical
//! (no transcendentals here).

/// Port of `phase/HmmUpdater.java` (static methods).
pub struct HmmUpdater;

impl HmmUpdater {
    /// `fwdUpdate(float[] fwd, float fwdSum, float pSwitch, float[] pMismatch, byte[] mismatch,
    /// int nStates)` — updates `fwd` in place and returns the new sum.
    pub fn fwd_update(
        fwd: &mut [f32],
        fwd_sum: f32,
        p_switch: f32,
        p_mismatch: &[f32],
        mismatch: &[u8],
        n_states: i32,
    ) -> f32 {
        assert!(p_mismatch.len() == 2, "{}", p_mismatch.len());
        let shift = p_switch / n_states as f32;
        let scale = (1.0f32 - p_switch) / fwd_sum;
        let mut new_sum = 0.0f32;
        for k in 0..n_states as usize {
            fwd[k] = p_mismatch[mismatch[k] as usize] * (scale * fwd[k] + shift);
            new_sum += fwd[k];
        }
        new_sum
    }

    /// `bwdUpdate(float[] bwd, float pSwitch, float[] pMismatch, byte[] mismatch, int nStates)`
    /// — updates `bwd` in place.
    pub fn bwd_update(
        bwd: &mut [f32],
        p_switch: f32,
        p_mismatch: &[f32],
        mismatch: &[u8],
        n_states: i32,
    ) {
        assert!(p_mismatch.len() == 2, "{}", p_mismatch.len());
        let mut sum = 0.0f32;
        for k in 0..n_states as usize {
            bwd[k] *= p_mismatch[mismatch[k] as usize];
            sum += bwd[k];
        }
        let shift = p_switch / n_states as f32;
        let scale = (1.0f32 - p_switch) / sum;
        for b in bwd.iter_mut().take(n_states as usize) {
            *b = scale * *b + shift;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fwd_update_matches_formula() {
        let mut fwd = [1.0f32, 1.0f32];
        let p_mismatch = [0.99f32, 0.01f32];
        let mismatch = [0u8, 1u8];
        let sum = HmmUpdater::fwd_update(&mut fwd, 2.0, 0.1, &p_mismatch, &mismatch, 2);
        // shift=0.05, scale=0.45; fwd0=0.99*0.5=0.495, fwd1=0.01*0.5=0.005
        assert!((fwd[0] - 0.495).abs() < 1e-6);
        assert!((fwd[1] - 0.005).abs() < 1e-6);
        assert!((sum - 0.5).abs() < 1e-6);
    }

    #[test]
    fn bwd_update_matches_formula() {
        let mut bwd = [1.0f32, 1.0f32];
        let p_mismatch = [0.99f32, 0.01f32];
        let mismatch = [0u8, 1u8];
        HmmUpdater::bwd_update(&mut bwd, 0.1, &p_mismatch, &mismatch, 2);
        // after *=pMismatch: [0.99, 0.01], sum=1.0; shift=0.05, scale=0.9
        // bwd0=0.9*0.99+0.05=0.941, bwd1=0.9*0.01+0.05=0.059
        assert!((bwd[0] - 0.941).abs() < 1e-6);
        assert!((bwd[1] - 0.059).abs() < 1e-6);
    }
}
