//! Port of the Java `phase` package — the haplotype phasing engine (PBWT-based IBS,
//! HMM forward/backward, and parameter estimation). Ported bottom-up from leaf types.

mod sample_seg;

pub use sample_seg::SampleSeg;
