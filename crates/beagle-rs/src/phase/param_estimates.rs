//! Port of `phase/ParamEstimates.java` — estimates the allele-mismatch probability and
//! recombination intensity for the Li–Stephens HMM.
//!
//! Java accumulates data in thread-safe queues and **sorts** before summing so the `double`
//! sums are order-independent (repeatable). The Rust port keeps that sort-then-sum to match
//! the floating-point result; `RefCell` provides the shared-mutability the Java queues gave.

use std::cell::RefCell;

#[derive(Clone, Copy)]
struct RecombData {
    gen_distance: f64,
    switch_prob: f64,
}

#[derive(Clone, Copy)]
struct MismatchData {
    marker_cnt: i32,
    p_mismatch_sum: f64,
}

/// Port of `phase/ParamEstimates.java`.
pub struct ParamEstimates {
    switch_data: RefCell<Vec<RecombData>>,
    mismatch_data: RefCell<Vec<MismatchData>>,
}

impl Default for ParamEstimates {
    fn default() -> Self {
        Self::new()
    }
}

impl ParamEstimates {
    /// `new ParamEstimates()`.
    pub fn new() -> Self {
        ParamEstimates {
            switch_data: RefCell::new(Vec::new()),
            mismatch_data: RefCell::new(Vec::new()),
        }
    }

    /// `addMismatchData(int markerCnt, double pMismatchSum)`.
    pub fn add_mismatch_data(&self, marker_cnt: i32, p_mismatch_sum: f64) {
        if marker_cnt > 0 && p_mismatch_sum > 0.0 && p_mismatch_sum.is_finite() {
            self.mismatch_data.borrow_mut().push(MismatchData {
                marker_cnt,
                p_mismatch_sum,
            });
        }
    }

    /// `addSwitchData(double genDistances, double switchProbs)`.
    pub fn add_switch_data(&self, gen_distances: f64, switch_probs: f64) {
        if gen_distances > 0.0
            && switch_probs > 0.0
            && gen_distances.is_finite()
            && switch_probs.is_finite()
        {
            self.switch_data.borrow_mut().push(RecombData {
                gen_distance: gen_distances,
                switch_prob: switch_probs,
            });
        }
    }

    /// `pMismatch()` — the estimated allele-mismatch rate, or `NaN` if no data.
    pub fn p_mismatch(&self) -> f32 {
        let mut mda = self.mismatch_data.borrow().clone();
        // compareTo: by pMismatchSum, then markerCnt
        mda.sort_by(|a, b| {
            a.p_mismatch_sum
                .total_cmp(&b.p_mismatch_sum)
                .then(a.marker_cnt.cmp(&b.marker_cnt))
        });
        let mut sum_markers: i64 = 0;
        let mut sum_p_mismatch = 0.0f64;
        for md in &mda {
            sum_markers += md.marker_cnt as i64;
            sum_p_mismatch += md.p_mismatch_sum;
        }
        if sum_markers == 0 {
            f32::NAN
        } else {
            (sum_p_mismatch / sum_markers as f64) as f32
        }
    }

    /// `recombIntensity()` — the estimated recombination intensity, or `NaN` if no data.
    pub fn recomb_intensity(&self) -> f32 {
        let mut rda = self.switch_data.borrow().clone();
        // compareTo: by genDistance, then switchProb
        rda.sort_by(|a, b| {
            a.gen_distance
                .total_cmp(&b.gen_distance)
                .then(a.switch_prob.total_cmp(&b.switch_prob))
        });
        let mut sum_switches = 0.0f64;
        let mut sum_distances = 0.0f64;
        for rd in &rda {
            sum_switches += rd.switch_prob;
            sum_distances += rd.gen_distance;
        }
        if sum_distances == 0.0 {
            f32::NAN
        } else {
            (sum_switches / sum_distances) as f32
        }
    }

    /// `clear()`.
    pub fn clear(&self) {
        self.switch_data.borrow_mut().clear();
        self.mismatch_data.borrow_mut().clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mismatch_rate_is_weighted_mean() {
        let pe = ParamEstimates::new();
        assert!(pe.p_mismatch().is_nan()); // no data
        pe.add_mismatch_data(100, 1.0);
        pe.add_mismatch_data(300, 3.0);
        // (1.0 + 3.0) / (100 + 300) = 0.01
        assert!((pe.p_mismatch() - 0.01).abs() < 1e-7);
        // non-positive/non-finite data is ignored
        pe.add_mismatch_data(0, 5.0);
        pe.add_mismatch_data(10, -1.0);
        assert!((pe.p_mismatch() - 0.01).abs() < 1e-7);
    }

    #[test]
    fn recomb_intensity_is_ratio() {
        let pe = ParamEstimates::new();
        assert!(pe.recomb_intensity().is_nan());
        pe.add_switch_data(2.0, 0.5);
        pe.add_switch_data(3.0, 1.0);
        // (0.5 + 1.0) / (2.0 + 3.0) = 0.3
        assert!((pe.recomb_intensity() - 0.3).abs() < 1e-7);
        pe.clear();
        assert!(pe.recomb_intensity().is_nan());
    }
}
