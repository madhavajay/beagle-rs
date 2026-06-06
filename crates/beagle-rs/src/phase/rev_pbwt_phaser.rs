//! Port of `phase/RevPbwtPhaser.java` — phases and imputes genotypes with the PBWT
//! processing markers in *decreasing* index order; heterozygotes the PBWT can't phase are
//! randomly phased and un-imputable alleles randomly imputed from the allele-frequency CDF.

use std::rc::Rc;

use crate::blbutil::BitArray;
use crate::jdk::Random;
use crate::vcf::GT;

use super::{FixedPhaseData, PbwtRecPhaser};

/// Port of `phase/RevPbwtPhaser.java`.
pub struct RevPbwtPhaser {
    targ_gt: Rc<dyn GT>,
    start: i32,
    end: i32,
    bits_per_allele: Vec<i32>,
    marker_to_bits: Vec<BitArray>,
}

fn impute_allele(allele_cdf: &[i32], rand: &mut Random) -> i32 {
    let bound = allele_cdf[allele_cdf.len() - 1];
    if bound == 0 {
        0
    } else {
        let r = rand.next_int_bound(allele_cdf[allele_cdf.len() - 1]);
        let mut allele = 0;
        while r >= allele_cdf[allele as usize] {
            allele += 1;
        }
        allele
    }
}

fn finish_phasing(
    alleles: &mut [i32],
    unph_het: &mut [bool],
    allele_cdf: &[i32],
    rand: &mut Random,
) {
    #[allow(clippy::needless_range_loop)] // s indexes unph_het and derives hap indices
    for s in 0..unph_het.len() {
        let h1 = s << 1;
        let h2 = h1 | 0b1;
        if unph_het[s] {
            let a1 = alleles[h1];
            let a2 = alleles[h2];
            if rand.next_boolean() {
                alleles[h1] = a2;
                alleles[h2] = a1;
            }
            unph_het[s] = false;
        } else {
            if alleles[h1] == -1 {
                alleles[h1] = impute_allele(allele_cdf, rand);
            }
            if alleles[h2] == -1 {
                alleles[h2] = impute_allele(allele_cdf, rand);
            }
        }
    }
}

fn store_phasing(targ_gt: &dyn GT, m: i32, alleles: &[i32]) -> BitArray {
    let n_targ_haps = targ_gt.n_haps();
    let bits_per_allele = targ_gt.markers().marker(m).bits_per_allele();
    let mut bits = BitArray::new(n_targ_haps * bits_per_allele);
    let mut bit = 0;
    #[allow(clippy::needless_range_loop)] // h is the haplotype index into alleles
    for h in 0..n_targ_haps as usize {
        let mut mask = 1;
        for _ in 0..bits_per_allele {
            if (alleles[h] & mask) == mask {
                bits.set(bit);
            }
            bit += 1;
            mask <<= 1;
        }
    }
    bits
}

fn phase(fpd: &FixedPhaseData, start: i32, end: i32, seed: i64) -> Vec<BitArray> {
    let mut rand = Random::new(seed);
    let overlap = fpd.stage1_overlap();
    let targ_gt = fpd.stage1_targ_gt().clone();
    let mut rec_phaser = PbwtRecPhaser::new(fpd);
    let mut missing_gt = vec![false; targ_gt.n_samples() as usize];
    let mut unph_het = vec![false; targ_gt.n_samples() as usize];
    let mut alleles = vec![0i32; fpd.n_haps() as usize];

    let mut marker_to_bits: Vec<Option<BitArray>> = (0..(end - start)).map(|_| None).collect();
    let mut last_m = -1;
    let mut m = end - 1;
    while m >= start {
        let allele_cdf = rec_phaser.phase(last_m, &mut alleles, m, &mut missing_gt, &mut unph_het);
        if m >= overlap {
            finish_phasing(&mut alleles, &mut unph_het, &allele_cdf, &mut rand);
        }
        marker_to_bits[(m - start) as usize] = Some(store_phasing(targ_gt.as_ref(), m, &alleles));
        last_m = m;
        m -= 1;
    }
    marker_to_bits
        .into_iter()
        .map(|b| b.expect("phased"))
        .collect()
}

impl RevPbwtPhaser {
    /// `new RevPbwtPhaser(FixedPhaseData fpd, int start, int end, long seed)`.
    pub fn new(fpd: &FixedPhaseData, start: i32, end: i32, seed: i64) -> Self {
        assert!(
            start >= 0 && end <= fpd.stage1_targ_gt().n_markers() && start < end,
            "{start}"
        );
        let targ_gt = fpd.stage1_targ_gt().clone();
        let bits_per_allele: Vec<i32> = (start..end)
            .map(|m| targ_gt.markers().marker(m).bits_per_allele())
            .collect();
        let marker_to_bits = phase(fpd, start, end, seed);
        RevPbwtPhaser {
            targ_gt,
            start,
            end,
            bits_per_allele,
            marker_to_bits,
        }
    }

    /// `targGT()`.
    pub fn targ_gt(&self) -> &Rc<dyn GT> {
        &self.targ_gt
    }
    /// `start()`.
    pub fn start(&self) -> i32 {
        self.start
    }
    /// `end()`.
    pub fn end(&self) -> i32 {
        self.end
    }
    /// `bitsPerAllele(int marker)`.
    pub fn bits_per_allele(&self, marker: i32) -> i32 {
        self.bits_per_allele[(marker - self.start) as usize]
    }

    /// `allele(int marker, int hap)`.
    pub fn allele(&self, marker: i32, hap: i32) -> i32 {
        let n_bits_per_allele = self.bits_per_allele[(marker - self.start) as usize];
        let b_start = hap * n_bits_per_allele;
        let bits = &self.marker_to_bits[(marker - self.start) as usize];
        if n_bits_per_allele == 1 {
            return i32::from(bits.get(b_start));
        }
        let mut allele = 0;
        let mut mask = 1;
        let b_end = b_start + n_bits_per_allele;
        for j in b_start..b_end {
            if bits.get(j) {
                allele |= mask;
            }
            mask <<= 1;
        }
        allele
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, HEADER_PREFIX,
    };
    use std::io::Write;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_rev_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..3 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 3],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let recs: Vec<Rc<dyn GTRec>> = (0..4)
            .map(|i| {
                let pos = 1_000_000 + i * 1_000_000;
                let line = format!("chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1");
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(
                    &h, &line, &mp,
                ))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(recs);
        let n = targ.n_markers();
        let indices = MarkerIndices::from_counts(0, n, n);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        let w = Window::new(gen_map, 1, true, indices, None, targ);
        let ped = Pedigree::new(w.targ_gt().samples().clone(), None);
        FixedPhaseData::new(&par, &ped, &w, None)
    }

    #[test]
    fn phases_reverse_and_preserves_genotypes() {
        let fpd = fpd();
        let n_markers = fpd.stage1_targ_gt().n_markers();
        let rp = RevPbwtPhaser::new(&fpd, 0, n_markers, 99999);
        assert_eq!(rp.start(), 0);
        assert_eq!(rp.end(), n_markers);
        // every marker's stored alleles form a valid het/hom matching the unordered input:
        // S0=0|0, S1=0/1 (some order), S2=1|1 at each marker
        for m in 0..n_markers {
            assert_eq!(rp.allele(m, 0), 0); // S0 hap0
            assert_eq!(rp.allele(m, 1), 0); // S0 hap1
            let s1 = (rp.allele(m, 2), rp.allele(m, 3));
            assert!(s1 == (0, 1) || s1 == (1, 0)); // S1 het, some phase
            assert_eq!(rp.allele(m, 4), 1); // S2 hap0
            assert_eq!(rp.allele(m, 5), 1); // S2 hap1
        }
    }
}
