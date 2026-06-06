//! Port of `vcf/BasicGTRec.java` — genotypes for one marker stored as a plain allele
//! array plus per-sample phased flags.
//!
//! The `VcfRecGTParser`-based constructor is added when `VcfRecGTParser` is ported; the
//! direct constructor here builds the same state from its components.

use crate::ints::IntArray;

use super::{to_vcf_rec, GTRec, Marker, Samples, VcfRecGTParser};

/// Port of `vcf/BasicGTRec.java`.
pub struct BasicGTRec {
    marker: Marker,
    samples: Samples,
    alleles: Vec<i32>,
    is_phased: Vec<bool>,
    all_phased: bool,
}

impl BasicGTRec {
    /// Builds a `BasicGTRec` from its components (`alleles` has length `2*nSamples`,
    /// `is_phased`/`all_phased` mirror `VcfRecGTParser.storeAlleles`).
    pub fn new(
        marker: Marker,
        samples: Samples,
        alleles: Vec<i32>,
        is_phased: Vec<bool>,
        all_phased: bool,
    ) -> Self {
        debug_assert_eq!(alleles.len(), 2 * samples.size() as usize);
        debug_assert_eq!(is_phased.len(), samples.size() as usize);
        BasicGTRec {
            marker,
            samples,
            alleles,
            is_phased,
            all_phased,
        }
    }

    /// `new BasicGTRec(VcfRecGTParser recParser)`.
    pub fn from_parser(rec_parser: &VcfRecGTParser) -> Self {
        let n_samples = rec_parser.samples().size();
        let mut alleles = vec![0i32; (n_samples << 1) as usize];
        let mut is_phased = vec![false; n_samples as usize];
        let all_phased = rec_parser.store_alleles_int(&mut alleles, &mut is_phased);
        BasicGTRec {
            marker: rec_parser.marker().clone(),
            samples: rec_parser.samples().clone(),
            alleles,
            is_phased,
            all_phased,
        }
    }
}

impl IntArray for BasicGTRec {
    fn size(&self) -> i32 {
        2 * self.samples.size()
    }

    fn get(&self, hap: i32) -> i32 {
        self.alleles[hap as usize]
    }
}

impl GTRec for BasicGTRec {
    fn samples(&self) -> &Samples {
        &self.samples
    }

    fn marker(&self) -> &Marker {
        &self.marker
    }

    fn is_phased_sample(&self, sample: i32) -> bool {
        self.is_phased[sample as usize]
    }

    fn is_phased(&self) -> bool {
        self.all_phased
    }
}

impl std::fmt::Display for BasicGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    #[test]
    fn basic_record() {
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &MarkerParser::new(true, true, true, true),
        );
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        let rec = BasicGTRec::new(marker, samples, vec![0, 1, 1, 0], vec![true, false], false);
        assert_eq!(rec.size(), 4);
        assert_eq!(rec.get(0), 0);
        assert_eq!(rec.get(3), 0);
        assert!(rec.is_phased_sample(0));
        assert!(!rec.is_phased_sample(1));
        assert!(!rec.is_phased());
        assert_eq!(
            rec.to_string(),
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1/0"
        );
    }
}
