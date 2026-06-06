//! Port of the Java `phase` package — the haplotype phasing engine (PBWT-based IBS,
//! HMM forward/backward, and parameter estimation). Ported bottom-up from leaf types.

mod fixed_phase_data;
mod hmm_updater;
mod ibs2;
mod ibs2_markers;
mod ibs2_sets;
mod param_estimates;
mod sample_phase;
mod sample_seg;
mod swap_rate;

pub use fixed_phase_data::FixedPhaseData;
pub use hmm_updater::HmmUpdater;
pub use ibs2::Ibs2;
pub use ibs2_markers::Ibs2Markers;
pub use ibs2_sets::Ibs2Sets;
pub use param_estimates::ParamEstimates;
pub use sample_phase::{ClustType, SamplePhase};
pub use sample_seg::SampleSeg;
pub use swap_rate::SwapRate;
