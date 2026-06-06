//! Port of `phase/FwdPbwtPhaser.java` — PBWT phasing in *increasing* marker order;
//! heterozygotes/missing alleles the forward PBWT can't resolve are reconciled against the
//! reverse PBWT phaser (and ultimately its random fallback).

use std::rc::Rc;

use crate::blbutil::BitArray;
use crate::vcf::GT;

use super::{FixedPhaseData, PbwtRecPhaser, RevPbwtPhaser};

/// Port of `phase/FwdPbwtPhaser.java`.
pub struct FwdPbwtPhaser {
    targ_gt: Rc<dyn GT>,
    start: i32,
    end: i32,
    bits_per_allele: Vec<i32>,
    marker_to_bits: Vec<BitArray>,
}

fn unpack_allele(bits: &BitArray, hap: i32, bits_per_allele: i32) -> i32 {
    let b_start = hap * bits_per_allele;
    if bits_per_allele == 1 {
        return i32::from(bits.get(b_start));
    }
    let b_end = b_start + bits_per_allele;
    let mut allele = 0;
    let mut mask = 1;
    for j in b_start..b_end {
        if bits.get(j) {
            allele |= mask;
        }
        mask <<= 1;
    }
    allele
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

fn update_last_het(alleles: &[i32], missing_gt: &[bool], last_het: &mut [i32], m: i32) {
    #[allow(clippy::needless_range_loop)] // s indexes missing_gt and derives hap indices
    for s in 0..missing_gt.len() {
        let h1 = s << 1;
        let h2 = h1 | 0b1;
        if !missing_gt[s] && alleles[h1] != alleles[h2] {
            last_het[s] = m;
        }
    }
}

fn impute_allele(
    marker_to_bits: &[Option<BitArray>],
    rev_pbwt: &RevPbwtPhaser,
    start: i32,
    last_het: i32,
    m: i32,
    hap: i32,
) -> i32 {
    if last_het < 0 {
        return rev_pbwt.allele(m, hap);
    }
    let comp_hap = hap ^ 0b1;
    let a1 = rev_pbwt.allele(last_het, hap);
    let a2 = rev_pbwt.allele(last_het, comp_hap);
    let bits_per_allele = rev_pbwt.bits_per_allele(last_het);
    let bits = marker_to_bits[(last_het - start) as usize]
        .as_ref()
        .expect("stored");
    let b1 = unpack_allele(bits, hap, bits_per_allele);
    let b2 = unpack_allele(bits, comp_hap, bits_per_allele);
    if (a1 < a2) == (b1 < b2) {
        rev_pbwt.allele(m, hap)
    } else {
        rev_pbwt.allele(m, comp_hap)
    }
}

#[allow(clippy::too_many_arguments)]
fn finish_phasing(
    marker_to_bits: &[Option<BitArray>],
    rev_pbwt: &RevPbwtPhaser,
    start: i32,
    m: i32,
    alleles: &mut [i32],
    last_het: &[i32],
    unph_het: &mut [bool],
) {
    #[allow(clippy::needless_range_loop)] // s indexes unph_het and derives hap indices
    for s in 0..unph_het.len() {
        let h1 = s << 1;
        let h2 = h1 | 0b1;
        if unph_het[s] {
            let prev_het = last_het[s];
            if prev_het >= 0 {
                let a1 = rev_pbwt.allele(prev_het, h1 as i32);
                let a2 = rev_pbwt.allele(prev_het, h2 as i32);
                let b1 = rev_pbwt.allele(m, h1 as i32);
                let b2 = rev_pbwt.allele(m, h2 as i32);
                let rev_same_phase = (a1 < a2) == (b1 < b2);
                let bits_per_allele = rev_pbwt.bits_per_allele(prev_het);
                let bits = marker_to_bits[(prev_het - start) as usize]
                    .as_ref()
                    .expect("stored");
                let c1 = unpack_allele(bits, h1 as i32, bits_per_allele);
                let c2 = unpack_allele(bits, h2 as i32, bits_per_allele);
                let fwd_same_phase = (c1 < c2) == (alleles[h1] < alleles[h2]);
                if rev_same_phase != fwd_same_phase {
                    alleles.swap(h1, h2);
                }
            }
            unph_het[s] = false;
        } else {
            if alleles[h1] == -1 {
                alleles[h1] =
                    impute_allele(marker_to_bits, rev_pbwt, start, last_het[s], m, h1 as i32);
            }
            if alleles[h2] == -1 {
                alleles[h2] =
                    impute_allele(marker_to_bits, rev_pbwt, start, last_het[s], m, h2 as i32);
            }
        }
    }
}

fn phase(fpd: &FixedPhaseData, start: i32, end: i32, seed: i64) -> Vec<BitArray> {
    let overlap = fpd.stage1_overlap();
    let targ_gt = fpd.stage1_targ_gt().clone();
    let mut rec_phaser = PbwtRecPhaser::new(fpd);
    let rev_pbwt = RevPbwtPhaser::new(fpd, start, end, seed);
    let n_samples = targ_gt.n_samples() as usize;
    let mut missing_gt = vec![false; n_samples];
    let mut unph_het = vec![false; n_samples];
    let mut last_het = vec![-1i32; n_samples];
    let mut alleles = vec![0i32; fpd.n_haps() as usize];

    let mut mkr_to_bits: Vec<Option<BitArray>> = (0..(end - start)).map(|_| None).collect();
    let mut last_m = -1;
    for m in start..end {
        rec_phaser.phase(last_m, &mut alleles, m, &mut missing_gt, &mut unph_het);
        if m >= overlap {
            finish_phasing(
                &mkr_to_bits,
                &rev_pbwt,
                start,
                m,
                &mut alleles,
                &last_het,
                &mut unph_het,
            );
        }
        mkr_to_bits[(m - start) as usize] = Some(store_phasing(targ_gt.as_ref(), m, &alleles));
        update_last_het(&alleles, &missing_gt, &mut last_het, m);
        last_m = m;
    }
    mkr_to_bits
        .into_iter()
        .map(|b| b.expect("phased"))
        .collect()
}

impl FwdPbwtPhaser {
    /// `new FwdPbwtPhaser(FixedPhaseData fpd, int start, int end, long seed)`.
    pub fn new(fpd: &FixedPhaseData, start: i32, end: i32, seed: i64) -> Self {
        assert!(
            start >= 0 && end <= fpd.targ_gt().n_markers() && start < end,
            "{start}"
        );
        let targ_gt = fpd.stage1_targ_gt().clone();
        let bits_per_allele: Vec<i32> = (start..end)
            .map(|m| targ_gt.markers().marker(m).bits_per_allele())
            .collect();
        let marker_to_bits = phase(fpd, start, end, seed);
        FwdPbwtPhaser {
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
    /// `allele(int marker, int hap)`.
    pub fn allele(&self, marker: i32, hap: i32) -> i32 {
        let index = (marker - self.start) as usize;
        unpack_allele(
            &self.marker_to_bits[index],
            hap,
            self.bits_per_allele[index],
        )
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
        let gt = std::env::temp_dir().join("beagle_rs_fwd_gt.vcf");
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
    fn forward_phasing_preserves_genotypes() {
        let fpd = fpd();
        let n = fpd.stage1_targ_gt().n_markers();
        let fp = FwdPbwtPhaser::new(&fpd, 0, n, 99999);
        assert_eq!((fp.start(), fp.end()), (0, n));
        for m in 0..n {
            assert_eq!(fp.allele(m, 0), 0);
            assert_eq!(fp.allele(m, 1), 0);
            let s1 = (fp.allele(m, 2), fp.allele(m, 3));
            assert!(s1 == (0, 1) || s1 == (1, 0));
            assert_eq!(fp.allele(m, 4), 1);
            assert_eq!(fp.allele(m, 5), 1);
        }
    }
}
