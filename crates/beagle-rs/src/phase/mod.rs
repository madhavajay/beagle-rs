//! Port of the Java `phase` package — the haplotype phasing engine (PBWT-based IBS,
//! HMM forward/backward, and parameter estimation). Ported bottom-up from leaf types.

mod hmm_updater;
mod param_estimates;
mod sample_phase;
mod sample_seg;
mod swap_rate;

pub use hmm_updater::HmmUpdater;
pub use param_estimates::ParamEstimates;
pub use sample_phase::{ClustType, SamplePhase};
pub use sample_seg::SampleSeg;
pub use swap_rate::SwapRate;
