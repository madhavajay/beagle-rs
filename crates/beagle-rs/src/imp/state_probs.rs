//! Port of `imp/StateProbs.java` — the (immutable) interface for a subset of Li & Stephens HMM
//! states and their probabilities for one target haplotype.

/// Port of the `imp/StateProbs.java` interface.
pub trait StateProbs {
    /// `targHap()`.
    fn targ_hap(&self) -> i32;

    /// `nTargMarkers()`.
    fn n_targ_markers(&self) -> i32;

    /// `nStates(int targMarker)` — number of stored states at `targMarker`.
    fn n_states(&self, targ_marker: i32) -> i32;

    /// `refHap(int targMarker, int index)`.
    fn ref_hap(&self, targ_marker: i32, index: i32) -> i32;

    /// `probs(int targMarker, int index)`.
    fn probs(&self, targ_marker: i32, index: i32) -> f32;

    /// `probsP1(int targMarker, int index)` — probability at the next marker (or this marker's
    /// probability when `targMarker + 1 == nTargMarkers()`).
    fn probs_p1(&self, targ_marker: i32, index: i32) -> f32;
}
