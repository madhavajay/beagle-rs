//! Port of `vcf/GTRec.java` — genotype data for one marker (an `IntArray` of alleles,
//! one per haplotype) plus the marker/samples and phasing info, with static helpers.

use crate::blbutil::consts;
use crate::ints::IntArray;

use super::{marker_utils, Marker, Samples};

/// Port of the `vcf/GTRec.java` interface (`extends IntArray`).
pub trait GTRec: IntArray {
    /// `samples()`.
    fn samples(&self) -> &Samples;

    /// `marker()`.
    fn marker(&self) -> &Marker;

    /// `isPhased(int sample)`.
    fn is_phased_sample(&self, sample: i32) -> bool;

    /// `isPhased()` — true iff every sample's genotype is phased and non-missing.
    fn is_phased(&self) -> bool;
}

/// `GTRec.alleleCounts(rec)` — count of each allele over all haplotypes.
pub fn allele_counts(rec: &dyn GTRec) -> Vec<i32> {
    let n_alleles = rec.marker().n_alleles();
    let mut cnts = vec![0i32; n_alleles as usize];
    for h in 0..rec.size() {
        let allele = rec.get(h);
        if allele >= 0 {
            cnts[allele as usize] += 1;
        }
    }
    cnts
}

/// `GTRec.alleleFreq(rec)`.
pub fn allele_freq(rec: &dyn GTRec) -> Vec<f64> {
    let cnts = allele_counts(rec);
    let sum: i32 = cnts.iter().sum();
    let mut freq = vec![0.0f64; cnts.len()];
    if sum > 0 {
        for (al, &c) in cnts.iter().enumerate() {
            freq[al] = c as f64 / sum as f64;
        }
    }
    freq
}

/// `GTRec.toVcfRec(gtRec)` — a VCF record line (missing QUAL/INFO are from the marker;
/// FILTER from the marker; GT format).
pub fn to_vcf_rec(rec: &dyn GTRec) -> String {
    let marker = rec.marker();
    let mut sb = String::with_capacity(100);
    marker_utils::append_first_7_fields(marker, &mut sb);
    sb.push(consts::TAB);
    sb.push_str(&marker.info());
    sb.push(consts::TAB);
    sb.push_str("GT");
    let n = rec.samples().size();
    for s in 0..n {
        let hap1 = s << 1;
        let a1 = rec.get(hap1);
        let a2 = rec.get(hap1 | 0b1);
        sb.push(consts::TAB);
        if a1 == -1 {
            sb.push(consts::MISSING_DATA_CHAR);
        } else {
            sb.push_str(&a1.to_string());
        }
        sb.push(if rec.is_phased_sample(s) {
            consts::PHASED_SEP
        } else {
            consts::UNPHASED_SEP
        });
        if a2 == -1 {
            sb.push(consts::MISSING_DATA_CHAR);
        } else {
            sb.push_str(&a2.to_string());
        }
    }
    sb
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    // A minimal in-memory GTRec for exercising the static helpers.
    struct TestRec {
        marker: Marker,
        samples: Samples,
        alleles: Vec<i32>,
        phased: Vec<bool>,
    }

    impl IntArray for TestRec {
        fn size(&self) -> i32 {
            self.alleles.len() as i32
        }
        fn get(&self, index: i32) -> i32 {
            self.alleles[index as usize]
        }
    }
    impl GTRec for TestRec {
        fn samples(&self) -> &Samples {
            &self.samples
        }
        fn marker(&self) -> &Marker {
            &self.marker
        }
        fn is_phased_sample(&self, sample: i32) -> bool {
            self.phased[sample as usize]
        }
        fn is_phased(&self) -> bool {
            self.phased.iter().all(|&p| p)
        }
    }

    fn rec() -> TestRec {
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &MarkerParser::new(true, true, true, true),
        );
        TestRec {
            marker,
            samples: Samples::new(&["s0".into(), "s1".into()], &[true, true]),
            alleles: vec![0, 1, 1, -1], // s0=0|1, s1=1|. (missing)
            phased: vec![true, false],
        }
    }

    #[test]
    fn counts_and_freq() {
        let r = rec();
        assert_eq!(allele_counts(&r), vec![1, 2]); // one 0, two 1s, one missing
        let f = allele_freq(&r);
        assert!((f[0] - 1.0 / 3.0).abs() < 1e-12);
        assert!((f[1] - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn to_vcf_rec_format() {
        let r = rec();
        // missing allele -> '.', phased s0 uses '|', unphased s1 uses '/'
        assert_eq!(
            to_vcf_rec(&r),
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1/."
        );
    }
}
