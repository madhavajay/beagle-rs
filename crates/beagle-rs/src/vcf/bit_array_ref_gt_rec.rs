//! Port of `vcf/BitArrayRefGTRec.java` — phased, non-missing genotypes for one marker,
//! with haplotype alleles packed `bitsPerAllele` bits each (LSB first) into a `BitArray`.

use crate::blbutil::BitArray;
use crate::ints::IntArray;
use crate::phase::{EstPhase, SamplePhase};

use super::{to_vcf_rec, GTRec, Marker, Samples};

/// Port of `vcf/BitArrayRefGTRec.java`.
pub struct BitArrayRefGTRec {
    bits_per_allele: i32,
    marker: Marker,
    samples: Samples,
    alleles: BitArray,
}

impl BitArrayRefGTRec {
    /// `BitArrayRefGTRec.toBitArrayRefGTRecs(EstPhase)` — converts the current per-sample
    /// (column-major) phase estimates into one row-major `BitArrayRefGTRec` per marker.
    /// Java builds the array in parallel; the order is preserved, so the Rust port builds it
    /// sequentially.
    pub fn to_bit_array_ref_gt_recs(est_phase: &EstPhase) -> Vec<BitArrayRefGTRec> {
        let fpd = est_phase.fpd();
        let gt = fpd.stage1_targ_gt();
        let markers = gt.markers().clone();
        let samples = gt.samples().clone();
        let bit_lists = SamplePhase::to_bit_lists(est_phase);
        bit_lists
            .into_iter()
            .enumerate()
            .map(|(m, bits)| {
                BitArrayRefGTRec::new(markers.marker(m as i32).clone(), samples.clone(), bits)
            })
            .collect()
    }

    /// `new BitArrayRefGTRec(Marker, Samples, BitArray)`.
    pub fn new(marker: Marker, samples: Samples, alleles: BitArray) -> Self {
        BitArrayRefGTRec {
            bits_per_allele: marker.bits_per_allele(),
            marker,
            samples,
            alleles,
        }
    }

    fn allele(&self, hap: i32) -> i32 {
        let start = self.bits_per_allele * hap;
        let end = start + self.bits_per_allele;
        let mut allele = 0i32;
        let mut mask = 1i32;
        for j in start..end {
            if self.alleles.get(j) {
                allele += mask;
            }
            mask <<= 1;
        }
        allele
    }
}

impl IntArray for BitArrayRefGTRec {
    fn size(&self) -> i32 {
        self.samples.size() << 1
    }

    fn get(&self, hap: i32) -> i32 {
        self.allele(hap)
    }
}

impl GTRec for BitArrayRefGTRec {
    fn samples(&self) -> &Samples {
        &self.samples
    }
    fn marker(&self) -> &Marker {
        &self.marker
    }
    fn is_phased_sample(&self, _sample: i32) -> bool {
        true
    }
    fn is_phased(&self) -> bool {
        true
    }
}

impl std::fmt::Display for BitArrayRefGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    fn mk(rec: &str) -> Marker {
        Marker::instance(rec, &MarkerParser::new(true, true, true, true))
    }

    #[test]
    fn biallelic_one_bit_per_allele() {
        let marker = mk("chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1"); // 1 bit/allele
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        let alleles = [0, 1, 1, 0];
        let mut bits = BitArray::new(4); // 1 bit * 4 haps
        for (h, &a) in alleles.iter().enumerate() {
            if a == 1 {
                bits.set(h as i32);
            }
        }
        let rec = BitArrayRefGTRec::new(marker, samples, bits);
        assert_eq!(rec.size(), 4);
        assert_eq!(
            (0..4).map(|h| rec.get(h)).collect::<Vec<_>>(),
            vec![0, 1, 1, 0]
        );
        assert!(rec.is_phased());
    }

    #[test]
    fn triallelic_two_bits_per_allele() {
        let marker = mk("chr1\t100\t.\tA\tC,G\t.\tPASS\t.\tGT\t0|1"); // 2 bits/allele
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        let alleles = [0, 1, 2, 0]; // 4 haps, 2 bits each -> 8 bits
        let mut bits = BitArray::new(8);
        for (h, &a) in alleles.iter().enumerate() {
            let start = 2 * h as i32;
            if a & 1 == 1 {
                bits.set(start);
            }
            if a & 2 == 2 {
                bits.set(start + 1);
            }
        }
        let rec = BitArrayRefGTRec::new(marker, samples, bits);
        assert_eq!(
            (0..4).map(|h| rec.get(h)).collect::<Vec<_>>(),
            vec![0, 1, 2, 0]
        );
    }
}
