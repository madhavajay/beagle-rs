//! Port of `vcf/GT.java` — genotype data for a list of markers and samples.

use super::{Marker, Markers, Samples};

/// Port of the `vcf/GT.java` interface.
pub trait GT {
    /// `isReversed()` — markers in decreasing base-position order.
    fn is_reversed(&self) -> bool;

    /// `nMarkers()`.
    fn n_markers(&self) -> i32;

    /// `marker(int)`.
    fn marker(&self, marker: i32) -> &Marker;

    /// `markers()` — markers in increasing chromosome-position order.
    fn markers(&self) -> &Markers;

    /// `nHaps()` = `2 * nSamples()`.
    fn n_haps(&self) -> i32;

    /// `nSamples()`.
    fn n_samples(&self) -> i32;

    /// `samples()`.
    fn samples(&self) -> &Samples;

    /// `isPhased()`.
    fn is_phased(&self) -> bool;

    /// `allele(int marker, int hap)` — the allele, or -1 if missing.
    fn allele(&self, marker: i32, hap: i32) -> i32;

    /// `restrict(Markers markers, int[] indices)`.
    fn restrict(&self, markers: &Markers, indices: &[i32]) -> Box<dyn GT>;

    /// `restrict(int start, int end)`.
    fn restrict_range(&self, start: i32, end: i32) -> Box<dyn GT>;
}
