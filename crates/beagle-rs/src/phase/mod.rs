//! Port of the Java `phase` package — the haplotype phasing engine (PBWT-based IBS,
//! HMM forward/backward, and parameter estimation). Ported bottom-up from leaf types.

mod coded_steps;
mod est_phase;
mod fixed_phase_data;
mod fwd_pbwt_phaser;
mod hmm_updater;
mod ibs2;
mod ibs2_markers;
mod ibs2_sets;
mod param_estimates;
mod pbwt_phaser;
mod pbwt_rec_phaser;
mod phase_data;
mod rev_pbwt_phaser;
mod sample_phase;
mod sample_seg;
mod swap_rate;

pub use coded_steps::CodedSteps;
pub use est_phase::EstPhase;
pub use fixed_phase_data::FixedPhaseData;
pub use fwd_pbwt_phaser::FwdPbwtPhaser;
pub use hmm_updater::HmmUpdater;
pub use ibs2::Ibs2;
pub use ibs2_markers::Ibs2Markers;
pub use ibs2_sets::Ibs2Sets;
pub use param_estimates::ParamEstimates;
pub use pbwt_phaser::PbwtPhaser;
pub use pbwt_rec_phaser::PbwtRecPhaser;
pub use phase_data::PhaseData;
pub use rev_pbwt_phaser::RevPbwtPhaser;
pub use sample_phase::{ClustType, SamplePhase};
pub use sample_seg::SampleSeg;
pub use swap_rate::SwapRate;
