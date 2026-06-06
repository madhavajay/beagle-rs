//! Port of `phase/PhaseBaum.java` — the interface for an HMM that updates the estimated
//! genotype phase of specified samples.

/// Port of the `phase/PhaseBaum.java` interface.
pub trait PhaseBaum {
    /// `nTargSamples()`.
    fn n_targ_samples(&self) -> i32;

    /// `phase(int sample)` — estimate and store the phased haplotypes for `sample`.
    fn phase(&mut self, sample: i32);
}
