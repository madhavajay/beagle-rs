//! Port of `vcf/BitArrayGTRec.java` — target (non-reference) genotypes for one marker,
//! storing haplotype alleles and a per-sample missing flag in `BitArray`s. All genotypes
//! are treated as unphased if any sample is unphased or missing.

use crate::blbutil::BitArray;
use crate::ints::IntArray;

use super::{to_vcf_rec, GTRec, HapListRep, Marker, Samples, VcfRecGTParser};

/// Port of `vcf/BitArrayGTRec.java`.
pub struct BitArrayGTRec {
    bits_per_allele: i32,
    marker: Marker,
    samples: Samples,
    is_phased: bool,
    is_missing: BitArray,
    alleles: BitArray,
}

fn store_allele(alleles: &mut BitArray, hap: i32, bits_per_allele: i32, allele: i32) {
    let base = hap * bits_per_allele;
    // Java loops a `mask = 1 << k`; `(allele & mask) == mask` <=> `(allele >> k) & 1 == 1`.
    for k in 0..bits_per_allele {
        if (allele >> k) & 1 == 1 {
            alleles.set(base + k);
        }
    }
}

impl BitArrayGTRec {
    /// `new BitArrayGTRec(VcfRecGTParser recParser)`.
    pub fn from_parser(rec_parser: &VcfRecGTParser) -> Self {
        let n_samples = rec_parser.samples().size();
        let n_haps = n_samples << 1;
        let bits_per_allele = rec_parser.marker().bits_per_allele();
        let marker = rec_parser.marker().clone();
        let samples = rec_parser.samples().clone();
        let mut allele_list = BitArray::new(n_haps * bits_per_allele);
        let mut is_missing_list = BitArray::new(n_samples);
        let is_phased = rec_parser.store_alleles_bits(&mut allele_list, &mut is_missing_list);
        BitArrayGTRec {
            bits_per_allele,
            marker,
            samples,
            is_phased,
            is_missing: is_missing_list,
            alleles: allele_list,
        }
    }

    /// `new BitArrayGTRec(VcfRecGTParser.HapListRep hlr)`.
    pub fn from_hap_list_rep(hlr: &HapListRep) -> Self {
        let n_samples = hlr.samples().size();
        let n_haps = n_samples << 1;
        let bits_per_allele = hlr.marker().bits_per_allele();
        let marker = hlr.marker().clone();
        let samples = hlr.samples().clone();
        let is_phased = hlr.is_phased();
        let mut is_missing = BitArray::new(n_samples);
        let mut alleles = BitArray::new(n_haps * bits_per_allele);
        let set_major_to_null = false;
        let hap_lists = hlr.hap_lists(set_major_to_null);
        let missing_samples = hlr.missing_samples();
        for (al, list) in hap_lists.iter().enumerate() {
            if let Some(list) = list {
                for &h in list {
                    store_allele(&mut alleles, h, bits_per_allele, al as i32);
                }
            }
        }
        for s in missing_samples {
            is_missing.set(s);
        }
        BitArrayGTRec {
            bits_per_allele,
            marker,
            samples,
            is_phased,
            is_missing,
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

impl IntArray for BitArrayGTRec {
    fn size(&self) -> i32 {
        2 * self.samples.size()
    }

    fn get(&self, hap: i32) -> i32 {
        if self.is_missing.get(hap >> 1) {
            -1
        } else {
            self.allele(hap)
        }
    }
}

impl GTRec for BitArrayGTRec {
    fn samples(&self) -> &Samples {
        &self.samples
    }
    fn marker(&self) -> &Marker {
        &self.marker
    }
    fn is_phased_sample(&self, _sample: i32) -> bool {
        self.is_phased
    }
    fn is_phased(&self) -> bool {
        self.is_phased
    }
}

impl std::fmt::Display for BitArrayGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{MarkerParser, VcfHeader, HEADER_PREFIX};

    fn header(n_dip: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n_dip {
            hdr.push_str(&format!("\tS{s}"));
        }
        let lines = vec!["##fileformat=VCFv4.2".to_string(), hdr];
        let is_dip = vec![true; n_dip];
        VcfHeader::new_accept_all("src", &lines, &is_dip)
    }

    // The `from_parser` path goes through `VcfRecGTParser::store_alleles_bits`, which
    // carries Beagle's latent quirk vcf-4: it stores `a1` for BOTH haplotypes of every
    // sample (never `a2`). This path is dead in real Beagle (only the HapListRep
    // constructor is reached, via VcfIt), so the bug is never observed; the port
    // preserves it and these tests pin the bug-preserving behavior.
    #[test]
    fn from_parser_preserves_a1_for_both_haps_quirk() {
        let h = header(2);
        let rec = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0";
        let p = VcfRecGTParser::new(&h, rec, &MarkerParser::new(true, true, true, true));
        let r = BitArrayGTRec::from_parser(&p);
        assert_eq!(r.size(), 4);
        // vcf-4: a1 stored for both haps -> [0,0] for 0|1, [1,1] for 1|0 (a2 dropped)
        assert_eq!(
            (0..4).map(|h| r.get(h)).collect::<Vec<_>>(),
            vec![0, 0, 1, 1]
        );
        assert!(r.is_phased());
    }

    #[test]
    fn from_parser_unphased_with_missing() {
        let h = header(2);
        // s1 is missing -> whole record unphased, get() returns -1 for missing sample's haps
        let rec = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t.|.";
        let p = VcfRecGTParser::new(&h, rec, &MarkerParser::new(true, true, true, true));
        let r = BitArrayGTRec::from_parser(&p);
        assert!(!r.is_phased());
        // vcf-4: a1 (=0) stored for both of s0's haps
        assert_eq!(r.get(0), 0);
        assert_eq!(r.get(1), 0);
        assert_eq!(r.get(2), -1);
        assert_eq!(r.get(3), -1);
    }

    #[test]
    fn from_hap_list_rep_matches_parser() {
        let h = header(3);
        let rec = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t1|1\t./.";
        let p = VcfRecGTParser::new(&h, rec, &MarkerParser::new(true, true, true, true));
        let r = BitArrayGTRec::from_hap_list_rep(&p.hap_list_rep());
        // s2 missing -> haps 4,5 are -1; s1 carries allele 1 on both haps
        assert_eq!(
            (0..6).map(|hp| r.get(hp)).collect::<Vec<_>>(),
            vec![0, 0, 1, 1, -1, -1]
        );
        assert!(!r.is_phased());
    }
}
