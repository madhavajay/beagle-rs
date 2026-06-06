//! Port of `vcf/SlidingWindow.java` — the interface for a sliding window of VCF records.
//!
//! The Java static `takeFromQ`/`addToQ` helpers coordinate the producer thread with the
//! consumer via a blocking queue; the Rust ports produce windows synchronously in
//! `next_window`, so those helpers are unnecessary.

use crate::main_pkg::Pedigree;

use super::{GeneticMap, Samples, Window};

/// Port of `vcf/SlidingWindow.java`.
pub trait SlidingWindow {
    /// `targSamples()`.
    fn targ_samples(&self) -> &Samples;

    /// `ped()`.
    fn ped(&self) -> &Pedigree;

    /// `genMap()`.
    fn gen_map(&self) -> &dyn GeneticMap;

    /// `cumTargMarkers()`.
    fn cum_targ_markers(&self) -> i32;

    /// `cumMarkers()`.
    fn cum_markers(&self) -> i32;

    /// `nextWindow()` — the next window, or `None` if there are no more.
    fn next_window(&mut self) -> Option<Window>;
}
