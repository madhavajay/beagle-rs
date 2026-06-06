//! Port of `vcf/VcfRecGTParser.java` — parses the GT format field of a VCF record into
//! allele arrays / haplotype lists. If one allele of a diploid genotype is missing, both
//! are set to missing.

use crate::blbutil::BitArray;

use super::marker_utils;
use super::{Marker, MarkerParser, Samples, VcfHeader};

const UNPHASED: u8 = b'/';
const PHASED: u8 = b'|';

/// Port of `vcf/VcfRecGTParser.java`. Borrows the header and the VCF record string.
pub struct VcfRecGTParser<'a> {
    vcf_header: &'a VcfHeader,
    vcf_rec: &'a str,
    marker: Marker,
    n_alleles: i32,
    n_samples: i32,
    ninth_tab_pos: i32,
}

fn index_of_tab_from(bytes: &[u8], from: i32) -> i32 {
    if from < 0 || from as usize > bytes.len() {
        return -1;
    }
    bytes[from as usize..]
        .iter()
        .position(|&b| b == b'\t')
        .map_or(-1, |o| from + o as i32)
}

/// `alEnd1` — exclusive end of the first allele (stops at separator/colon/tab).
fn al_end1(bytes: &[u8], start: usize) -> i32 {
    assert!(start != bytes.len(), "genotype is missing allele separator");
    let mut i = start;
    while i < bytes.len() {
        let c = bytes[i];
        if c == UNPHASED || c == PHASED || c == b'\t' || c == b':' {
            return i as i32;
        }
        i += 1;
    }
    i as i32
}

/// `alEnd2` — exclusive end of the second allele (stops at colon/tab).
fn al_end2(bytes: &[u8], start: i32) -> i32 {
    let mut i = start as usize;
    while i < bytes.len() {
        let c = bytes[i];
        if c == b':' || c == b'\t' {
            return i as i32;
        }
        i += 1;
    }
    i as i32
}

impl<'a> VcfRecGTParser<'a> {
    /// `new VcfRecGTParser(VcfHeader, String vcfRec, MarkerParser)`.
    pub fn new(vcf_header: &'a VcfHeader, vcf_rec: &'a str, marker_parser: &MarkerParser) -> Self {
        assert!(vcf_header.n_samples() != 0, "nSamples==0");
        let marker = Marker::instance(vcf_rec, marker_parser);
        let n_alleles = marker.n_alleles();
        let n_samples = vcf_header.n_samples();
        let ninth_tab_pos = marker_utils::ninth_tab_pos(vcf_rec);
        VcfRecGTParser {
            vcf_header,
            vcf_rec,
            marker,
            n_alleles,
            n_samples,
            ninth_tab_pos,
        }
    }

    /// `marker()`.
    pub fn marker(&self) -> &Marker {
        &self.marker
    }

    /// `nAlleles()`.
    pub fn n_alleles(&self) -> i32 {
        self.n_alleles
    }

    /// `samples()`.
    pub fn samples(&self) -> &Samples {
        self.vcf_header.samples()
    }

    /// `nSamples()`.
    pub fn n_samples(&self) -> i32 {
        self.n_samples
    }

    fn parse_allele(&self, start: i32, end: i32) -> i32 {
        let bytes = self.vcf_rec.as_bytes();
        assert!(
            start != end,
            "ERROR: Missing sample allele: {}",
            self.vcf_rec
        );
        let al = if start + 1 == end {
            let c = bytes[start as usize];
            if c == b'.' {
                return -1;
            }
            (c as i32) - (b'0' as i32)
        } else {
            self.vcf_rec[start as usize..end as usize]
                .parse::<i32>()
                .expect("parsable allele")
        };
        assert!(
            al >= 0 && al < self.n_alleles,
            "ERROR: Invalid allele at character {} in record \"{}\t...\"",
            start,
            self.marker
        );
        al
    }

    fn haploid_diploid_error(&self, sample: i32, is_diploid: bool) -> ! {
        panic!(
            "Sample {} has an inconsistent number of alleles. The first genotype is {}, \
             but the genotype at position {}:{} is {}",
            self.vcf_header.samples().id(sample),
            if self.samples().is_diploid(sample) {
                "diploid"
            } else {
                "haploid"
            },
            self.marker.chrom(),
            self.marker.pos(),
            if is_diploid { "diploid" } else { "haploid" }
        );
    }

    /// `storeAlleles(int[] alleles, boolean[] isPhased)`.
    pub fn store_alleles_int(&self, alleles: &mut [i32], is_phased: &mut [bool]) -> bool {
        assert!(
            alleles.len() == (self.n_samples << 1) as usize,
            "{}",
            alleles.len()
        );
        assert!(
            is_phased.len() == self.n_samples as usize,
            "{}",
            is_phased.len()
        );
        let bytes = self.vcf_rec.as_bytes();
        let mut pos = self.ninth_tab_pos;
        let mut unfilt: i32 = -1;
        let mut all_phased = true;
        for s in 0..self.n_samples {
            assert!(pos != -1, "VCF record field-count error");
            let next_unfiltered = self.vcf_header.unfiltered_sample_index(s);
            unfilt += 1;
            while unfilt < next_unfiltered {
                pos = index_of_tab_from(bytes, pos + 1);
                assert!(pos != -1, "VCF record field-count error");
                unfilt += 1;
            }
            let al_start = (pos + 1) as usize;
            let end1 = al_end1(bytes, al_start);
            assert!(al_start as i32 != end1, "missing data for sample {s}");
            let end2 = al_end2(bytes, end1);
            let is_diploid = end1 != end2;
            if is_diploid != self.samples().is_diploid(s) {
                self.haploid_diploid_error(s, is_diploid);
            }
            let h1 = (s << 1) as usize;
            let h2 = h1 | 0b1;
            let mut a1 = self.parse_allele(al_start as i32, end1);
            let mut a2 = if end1 == end2 {
                a1
            } else {
                self.parse_allele(end1 + 1, end2)
            };
            if (a1 == -1) ^ (a2 == -1) {
                a1 = -1;
                a2 = -1;
            }
            alleles[h1] = a1;
            alleles[h2] = a2;
            is_phased[s as usize] = !is_diploid || bytes[end1 as usize] == PHASED;
            all_phased &= is_phased[s as usize];
            pos = index_of_tab_from(bytes, end2);
        }
        all_phased
    }

    /// `storeAlleles(BitArray alleles, BitArray isMissing)`.
    pub fn store_alleles_bits(&self, alleles: &mut BitArray, is_missing: &mut BitArray) -> bool {
        let bits_per_allele = self.marker.bits_per_allele();
        assert!(is_missing.size() == self.n_samples, "{}", is_missing.size());
        assert!(
            alleles.size() == (self.n_samples << 1) * bits_per_allele,
            "{}",
            alleles.size()
        );
        let bytes = self.vcf_rec.as_bytes();
        let mut pos = self.ninth_tab_pos;
        let mut unfilt: i32 = -1;
        let mut is_phased = true;
        for s in 0..self.n_samples {
            assert!(pos != -1, "VCF record field-count error");
            let next_unfiltered = self.vcf_header.unfiltered_sample_index(s);
            unfilt += 1;
            while unfilt < next_unfiltered {
                pos = index_of_tab_from(bytes, pos + 1);
                assert!(pos != -1, "VCF record field-count error");
                unfilt += 1;
            }
            let al_start = (pos + 1) as usize;
            let end1 = al_end1(bytes, al_start);
            assert!(al_start as i32 != end1, "missing data for sample {s}");
            let end2 = al_end2(bytes, end1);
            let is_diploid = end1 != end2;
            if is_diploid != self.samples().is_diploid(s) {
                self.haploid_diploid_error(s, is_diploid);
            }
            let a1 = self.parse_allele(al_start as i32, end1);
            let a2 = if end1 == end2 {
                a1
            } else {
                self.parse_allele(end1 + 1, end2)
            };
            if is_diploid && bytes[end1 as usize] == UNPHASED {
                is_phased = false;
            }
            if a1 == -1 || a2 == -1 {
                is_phased = false;
                is_missing.set(s);
            } else {
                let h1 = s << 1;
                let h2 = h1 | 0b1;
                store_allele(alleles, h1, bits_per_allele, a1);
                store_allele(alleles, h2, bits_per_allele, a1);
            }
            pos = index_of_tab_from(bytes, end2);
        }
        is_phased
    }

    /// `phasedAlleles()` — requires all genotypes phased and non-missing.
    fn phased_alleles(&self) -> Vec<i32> {
        let bytes = self.vcf_rec.as_bytes();
        let mut alleles = vec![0i32; 2 * self.n_samples as usize];
        let mut pos = self.ninth_tab_pos;
        let mut unfilt: i32 = -1;
        let mut hap = 0usize;
        for s in 0..self.n_samples {
            assert!(pos != -1, "VCF record field-count error");
            let next_unfiltered = self.vcf_header.unfiltered_sample_index(s);
            unfilt += 1;
            while unfilt < next_unfiltered {
                pos = index_of_tab_from(bytes, pos + 1);
                assert!(pos != -1, "VCF record field-count error");
                unfilt += 1;
            }
            let al_start = (pos + 1) as usize;
            let end1 = al_end1(bytes, al_start);
            assert!(al_start as i32 != end1, "missing data for sample {s}");
            let end2 = al_end2(bytes, end1);
            let is_diploid = end1 != end2;
            let a1 = self.parse_allele(al_start as i32, end1);
            let a2 = if end1 == end2 {
                a1
            } else {
                self.parse_allele(end1 + 1, end2)
            };
            if is_diploid != self.samples().is_diploid(s) {
                self.haploid_diploid_error(s, is_diploid);
            }
            assert!(
                !((is_diploid && bytes[end1 as usize] != PHASED) || a1 == -1 || a2 == -1),
                "unphased or missing genotype for sample {s}"
            );
            alleles[hap] = a1;
            hap += 1;
            alleles[hap] = a2;
            hap += 1;
            pos = index_of_tab_from(bytes, end2);
        }
        alleles
    }

    /// `nonMajRefIndices()` — per allele, the haplotypes carrying a non-major allele
    /// (`None` for the major allele).
    pub fn non_maj_ref_indices(&self) -> Vec<Option<Vec<i32>>> {
        let alleles = self.phased_alleles();
        let mut al_cnts = vec![0i32; self.n_alleles as usize];
        for &a in &alleles {
            al_cnts[a as usize] += 1;
        }
        let mut maj_al = 0usize;
        for j in 1..self.n_alleles as usize {
            if al_cnts[j] > al_cnts[maj_al] {
                maj_al = j;
            }
        }
        let mut non_maj: Vec<Option<Vec<i32>>> = (0..self.n_alleles as usize)
            .map(|al| {
                if al == maj_al {
                    None
                } else {
                    Some(Vec::with_capacity(al_cnts[al] as usize))
                }
            })
            .collect();
        for (j, &al) in alleles.iter().enumerate() {
            if al as usize != maj_al {
                non_maj[al as usize].as_mut().unwrap().push(j as i32);
            }
        }
        non_maj
    }

    /// `hapListRep()`.
    pub fn hap_list_rep(&self) -> HapListRep {
        let bytes = self.vcf_rec.as_bytes();
        let mut hap_lists: Vec<Vec<i32>> = (0..self.n_alleles).map(|_| Vec::new()).collect();
        let mut miss_list = Vec::new();
        let mut is_phased = true;
        let mut tab_index = self.ninth_tab_pos;
        let mut unfilt: i32 = -1;
        for s in 0..self.n_samples {
            assert!(tab_index != -1, "VCF record field-count error");
            let next_unfiltered = self.vcf_header.unfiltered_sample_index(s);
            unfilt += 1;
            while unfilt < next_unfiltered {
                tab_index = index_of_tab_from(bytes, tab_index + 1);
                assert!(tab_index != -1, "VCF record field-count error");
                unfilt += 1;
            }
            let al_start = (tab_index + 1) as usize;
            let end1 = al_end1(bytes, al_start);
            assert!(al_start as i32 != end1, "missing data for sample {s}");
            let end2 = al_end2(bytes, end1);
            let is_diploid = end1 != end2;
            is_phased &= !is_diploid || bytes[end1 as usize] == PHASED;
            if is_diploid != self.samples().is_diploid(s) {
                self.haploid_diploid_error(s, is_diploid);
            }
            let h1 = s << 1;
            let h2 = h1 | 0b1;
            let a1 = self.parse_allele(al_start as i32, end1);
            let a2 = if end1 == end2 {
                a1
            } else {
                self.parse_allele(end1 + 1, end2)
            };
            if a1 < 0 || a2 < 0 {
                miss_list.push(s);
                is_phased = false;
            } else {
                hap_lists[a1 as usize].push(h1);
                hap_lists[a2 as usize].push(h2);
            }
            tab_index = index_of_tab_from(bytes, end2);
        }
        HapListRep::new(
            self.marker.clone(),
            self.samples().clone(),
            miss_list,
            hap_lists,
            is_phased,
        )
    }
}

fn store_allele(alleles: &mut BitArray, hap: i32, bits_per_allele: i32, allele: i32) {
    let base = hap * bits_per_allele;
    for k in 0..bits_per_allele {
        if (allele >> k) & 1 == 1 {
            alleles.set(base + k);
        }
    }
}

fn maj_allele(hap_lists: &[Vec<i32>]) -> i32 {
    let mut maj = 0usize;
    for j in 1..hap_lists.len() {
        if hap_lists[j].len() > hap_lists[maj].len() {
            maj = j;
        }
    }
    maj as i32
}

/// Port of `VcfRecGTParser.HapListRep`.
pub struct HapListRep {
    marker: Marker,
    samples: Samples,
    missing_samples: Vec<i32>,
    hap_lists: Vec<Vec<i32>>,
    is_phased: bool,
    major_allele: i32,
}

impl HapListRep {
    fn new(
        marker: Marker,
        samples: Samples,
        missing_samples: Vec<i32>,
        hap_lists: Vec<Vec<i32>>,
        is_phased: bool,
    ) -> Self {
        let major_allele = maj_allele(&hap_lists);
        HapListRep {
            marker,
            samples,
            missing_samples,
            hap_lists,
            is_phased,
            major_allele,
        }
    }

    /// `nonmajorAlleleCnt()`.
    pub fn nonmajor_allele_cnt(&self) -> i32 {
        let mut cnt = -(self.hap_lists[self.major_allele as usize].len() as i32);
        for list in &self.hap_lists {
            cnt += list.len() as i32;
        }
        cnt
    }

    /// `missingSamples()`.
    pub fn missing_samples(&self) -> Vec<i32> {
        self.missing_samples.clone()
    }

    /// `majorAllele()`.
    pub fn major_allele(&self) -> i32 {
        self.major_allele
    }

    /// `hapLists(boolean setMajorToNull)`.
    pub fn hap_lists(&self, set_major_to_null: bool) -> Vec<Option<Vec<i32>>> {
        self.hap_lists
            .iter()
            .enumerate()
            .map(|(j, list)| {
                if !set_major_to_null || j as i32 != self.major_allele {
                    Some(list.clone())
                } else {
                    None
                }
            })
            .collect()
    }

    /// `marker()`.
    pub fn marker(&self) -> &Marker {
        &self.marker
    }

    /// `samples()`.
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// `isPhased()`.
    pub fn is_phased(&self) -> bool {
        self.is_phased
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::HEADER_PREFIX;

    fn header(n_dip: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n_dip {
            hdr.push_str(&format!("\tS{s}"));
        }
        let lines = vec!["##fileformat=VCFv4.2".to_string(), hdr];
        let is_dip = vec![true; n_dip];
        VcfHeader::new_accept_all("src", &lines, &is_dip)
    }

    #[test]
    fn store_alleles_int_phased_unphased_missing() {
        let h = header(3);
        let rec = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1/1\t./.";
        let p = VcfRecGTParser::new(&h, rec, &MarkerParser::new(true, true, true, true));
        let mut alleles = vec![0; 6];
        let mut is_phased = vec![false; 3];
        let all_phased = p.store_alleles_int(&mut alleles, &mut is_phased);
        assert_eq!(alleles, vec![0, 1, 1, 1, -1, -1]);
        assert_eq!(is_phased, vec![true, false, false]);
        assert!(!all_phased);
    }

    #[test]
    fn non_maj_ref_indices_phased() {
        let h = header(4);
        // alleles per hap: s0=0|0, s1=0|1, s2=1|0, s3=0|0 -> [0,0,0,1,1,0,0,0]; major 0
        let rec = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t0|0";
        let p = VcfRecGTParser::new(&h, rec, &MarkerParser::new(true, true, true, true));
        let non_maj = p.non_maj_ref_indices();
        assert_eq!(non_maj[0], None); // major allele 0
        assert_eq!(non_maj[1], Some(vec![3, 4])); // haps carrying allele 1
    }

    #[test]
    fn hap_list_rep_with_missing() {
        let h = header(3);
        let rec = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t1|1\t./.";
        let p = VcfRecGTParser::new(&h, rec, &MarkerParser::new(true, true, true, true));
        let rep = p.hap_list_rep();
        assert_eq!(rep.missing_samples(), vec![2]);
        assert!(!rep.is_phased()); // missing -> not phased
        assert_eq!(rep.major_allele(), 0); // allele0 has 2 haps, allele1 has 2 -> first max = 0
        let lists = rep.hap_lists(true);
        assert_eq!(lists[0], None);
        assert_eq!(lists[1], Some(vec![2, 3]));
    }
}
