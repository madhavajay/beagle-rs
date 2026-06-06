//! Port of the Java `imp` package — the genotype imputation engine (Li & Stephens HMM over a
//! reference panel, imputing ungenotyped markers). Ported bottom-up from leaf types.

mod haplotype_coder;
mod ref_hap_hash;
mod state_probs;
mod state_probs_factory;

pub use haplotype_coder::HaplotypeCoder;
pub use ref_hap_hash::RefHapHash;
pub use state_probs::StateProbs;
pub use state_probs_factory::{BasicStateProbs, StateProbsFactory};
