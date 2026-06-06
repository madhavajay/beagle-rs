//! Port of `phase/SwapRate.java` — global counters for the proportion of unphased
//! heterozygotes whose phase relative to the previous heterozygote was reversed.

use std::sync::atomic::{AtomicI64, Ordering};

static N_SWAPS: AtomicI64 = AtomicI64::new(0);
static N_UNPH_HETS: AtomicI64 = AtomicI64::new(0);

/// Port of `phase/SwapRate.java` (static methods).
pub struct SwapRate;

impl SwapRate {
    /// `getAndResetSwapRate()` — `nSwaps / nUnphHets`, then resets both counters to 0.
    pub fn get_and_reset_swap_rate() -> f64 {
        let rate =
            N_SWAPS.load(Ordering::SeqCst) as f64 / N_UNPH_HETS.load(Ordering::SeqCst) as f64;
        N_SWAPS.store(0, Ordering::SeqCst);
        N_UNPH_HETS.store(0, Ordering::SeqCst);
        rate
    }

    /// `increment(int nUnphHets, int nSwaps)`.
    pub fn increment(n_unph_hets: i32, n_swaps: i32) {
        N_SWAPS.fetch_add(n_swaps as i64, Ordering::SeqCst);
        N_UNPH_HETS.fetch_add(n_unph_hets as i64, Ordering::SeqCst);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn increment_then_get_and_reset() {
        // reset any prior state first (global counters)
        SwapRate::get_and_reset_swap_rate();
        SwapRate::increment(10, 3);
        SwapRate::increment(10, 1);
        // 4 swaps / 20 unphased hets = 0.2
        assert!((SwapRate::get_and_reset_swap_rate() - 0.2).abs() < 1e-12);
        // counters reset -> 0/0 = NaN
        assert!(SwapRate::get_and_reset_swap_rate().is_nan());
    }
}
