//! Port of the Java `imp` package — the genotype imputation engine (Li & Stephens HMM over a
//! reference panel, imputing ungenotyped markers). Ported bottom-up from leaf types.

mod haplotype_coder;
mod state_probs;

pub use haplotype_coder::HaplotypeCoder;
pub use state_probs::StateProbs;
