//! Port of `phase/PbwtRecPhaser.java` — partial phasing and imputation of genotypes using
//! the Positional Burrows–Wheeler Transform (Durbin 2014; Delaneau et al. 2019).
//!
//! The static `bitSets(...)` helper is added with the Fwd/Rev PBWT phasers that consume it
//! (its element type depends on their usage).

use std::rc::Rc;

use crate::beagleutil::PbwtUpdater;
use crate::vcf::{RefGT, GT};

use super::FixedPhaseData;

/// Port of `phase/PbwtRecPhaser.java`.
pub struct PbwtRecPhaser {
    targ_gt: Rc<dyn GT>,
    opt_ref: Option<RefGT>,
    phased_overlap: i32,
    n_targ_haps: i32,
    n_targ_samples: i32,
    n_haps: i32,
    a: Vec<i32>,
    inv_a: Vec<i32>,
    pbwt: PbwtUpdater,
}

fn phase_cnt_static(adjacent_allele: i32, a1: i32, a2: i32) -> i32 {
    if adjacent_allele == a1 {
        1
    } else if adjacent_allele == a2 {
        -1
    } else {
        0
    }
}

impl PbwtRecPhaser {
    /// `new PbwtRecPhaser(FixedPhaseData fpd)`.
    pub fn new(fpd: &FixedPhaseData) -> Self {
        let targ_gt = fpd.stage1_targ_gt().clone();
        let opt_ref = fpd.stage1_ref_gt().cloned();
        let n_targ_haps = targ_gt.n_haps();
        let n_targ_samples = targ_gt.n_samples();
        let n_haps = targ_gt.n_haps() + opt_ref.as_ref().map_or(0, |r| r.n_haps());
        PbwtRecPhaser {
            targ_gt,
            opt_ref,
            phased_overlap: fpd.stage1_overlap(),
            n_targ_haps,
            n_targ_samples,
            n_haps,
            a: (0..n_haps).collect(),
            inv_a: (0..n_haps).collect(),
            pbwt: PbwtUpdater::new(n_haps),
        }
    }

    /// `nHaps()`.
    pub fn n_haps(&self) -> i32 {
        self.n_haps
    }

    /// `targGT()`.
    pub fn targ_gt(&self) -> &Rc<dyn GT> {
        &self.targ_gt
    }

    /// `phase(int currentMkr, int[] alleles, int nextMkr, boolean[] missing,
    /// boolean[] unphHet)` — returns the allele-count CDF at `nextMkr`.
    pub fn phase(
        &mut self,
        current_mkr: i32,
        alleles: &mut [i32],
        next_mkr: i32,
        missing: &mut [bool],
        unph_het: &mut [bool],
    ) -> Vec<i32> {
        self.check_arrays(alleles, missing, unph_het);
        if current_mkr != -1 {
            let n_alleles = self.targ_gt.marker(current_mkr).n_alleles();
            self.pbwt.update_alleles(alleles, n_alleles, &mut self.a);
        }
        let al_cnts = self.set_alleles(next_mkr, alleles, unph_het, missing);
        if next_mkr >= self.phased_overlap {
            self.phase_haps(alleles, unph_het);
        }
        al_cnts
    }

    fn check_arrays(&self, alleles: &[i32], missing: &[bool], unph_het: &[bool]) {
        assert!(alleles.len() as i32 == self.n_haps, "{}", alleles.len());
        assert!(
            missing.len() as i32 == self.n_targ_samples,
            "{}",
            missing.len()
        );
        assert!(
            unph_het.len() as i32 == self.n_targ_samples,
            "{}",
            unph_het.len()
        );
    }

    fn phase_haps(&mut self, alleles: &mut [i32], unph_het: &mut [bool]) {
        self.set_inv_a();
        let mut threshold = 2;
        let mut change_made = true;
        while threshold > 0 || change_made {
            change_made = false;
            for s in 0..self.n_targ_samples {
                if unph_het[s as usize] {
                    change_made |= self.phase_sample(s, threshold, alleles, unph_het);
                } else {
                    let h1 = (s << 1) as usize;
                    let h2 = h1 | 0b1;
                    if alleles[h1] == -1 {
                        let v = self.impute(alleles, unph_het, self.inv_a[h1]);
                        alleles[h1] = v;
                        change_made |= v >= 0;
                    }
                    if alleles[h2] == -1 {
                        let v = self.impute(alleles, unph_het, self.inv_a[h2]);
                        alleles[h2] = v;
                        change_made |= v >= 0;
                    }
                }
            }
            if !change_made {
                threshold -= 1;
            }
        }
    }

    fn phase_sample(
        &self,
        s: i32,
        threshold: i32,
        alleles: &mut [i32],
        unph_het: &mut [bool],
    ) -> bool {
        let h1 = (s << 1) as usize;
        let h2 = h1 | 0b1;
        let a1 = alleles[h1];
        let a2 = alleles[h2];
        debug_assert!(a1 >= 0 && a2 >= 0 && a1 != a2);
        let cnt1 = self.phase_cnt(alleles, unph_het, self.inv_a[h1], a1, a2);
        let cnt2 = self.phase_cnt(alleles, unph_het, self.inv_a[h2], a2, a1);
        let cnt = cnt1 + cnt2;
        if cnt >= threshold {
            unph_het[s as usize] = false;
            return true;
        }
        if cnt <= -threshold {
            alleles[h1] = a2;
            alleles[h2] = a1;
            unph_het[s as usize] = false;
            return true;
        }
        false
    }

    fn phase_cnt(&self, alleles: &[i32], unphased_het: &[bool], ai: i32, a1: i32, a2: i32) -> i32 {
        let mut phase_cnt = 0;
        if ai > 0 {
            let h = self.a[(ai - 1) as usize];
            let s = h >> 1;
            if s as usize >= unphased_het.len() || !unphased_het[s as usize] {
                phase_cnt += phase_cnt_static(alleles[h as usize], a1, a2);
            }
        }
        if (ai + 1) < alleles.len() as i32 {
            let h = self.a[(ai + 1) as usize];
            let s = h >> 1;
            if s as usize >= unphased_het.len() || !unphased_het[s as usize] {
                phase_cnt += phase_cnt_static(alleles[h as usize], a1, a2);
            }
        }
        phase_cnt
    }

    fn impute(&self, input_alleles: &[i32], unphased_het: &[bool], ai: i32) -> i32 {
        let mut prev = -1;
        let mut next = -1;
        if ai > 0 {
            let h = self.a[(ai - 1) as usize];
            let s = h >> 1;
            if s as usize >= unphased_het.len() || !unphased_het[s as usize] {
                prev = input_alleles[h as usize];
            }
        }
        if (ai + 1) < self.a.len() as i32 {
            let h = self.a[(ai + 1) as usize];
            let s = h >> 1;
            if s as usize >= unphased_het.len() || !unphased_het[s as usize] {
                next = input_alleles[h as usize];
            }
        }
        if prev >= 0 && (prev == next || next < 0) {
            prev
        } else if prev < 0 && next >= 0 {
            next
        } else {
            -1
        }
    }

    fn set_inv_a(&mut self) {
        for j in 0..self.a.len() {
            self.inv_a[self.a[j] as usize] = j as i32;
        }
    }

    fn set_alleles(
        &self,
        m: i32,
        input_alleles: &mut [i32],
        unph_het: &mut [bool],
        missing: &mut [bool],
    ) -> Vec<i32> {
        let mut al_cnts = vec![0i32; self.targ_gt.marker(m).n_alleles() as usize];
        self.set_targ_alleles(m, input_alleles, unph_het, missing, &mut al_cnts);
        if self.opt_ref.is_some() {
            self.set_ref_alleles(m, input_alleles, &mut al_cnts);
        }
        for j in 1..al_cnts.len() {
            al_cnts[j] += al_cnts[j - 1];
        }
        al_cnts
    }

    fn set_targ_alleles(
        &self,
        m: i32,
        input_alleles: &mut [i32],
        unph_het: &mut [bool],
        missing: &mut [bool],
        al_cnts: &mut [i32],
    ) {
        for s in 0..self.n_targ_samples {
            let h1 = (s << 1) as usize;
            let h2 = h1 | 0b1;
            let a1 = self.targ_gt.allele(m, h1 as i32);
            let a2 = self.targ_gt.allele(m, h2 as i32);
            input_alleles[h1] = a1;
            input_alleles[h2] = a2;
            unph_het[s as usize] = m >= self.phased_overlap && (a1 >= 0 && a2 >= 0 && a1 != a2);
            missing[s as usize] = a1 < 0 || a2 < 0;
            if a1 >= 0 {
                al_cnts[a1 as usize] += 1;
            }
            if a2 >= 0 {
                al_cnts[a2 as usize] += 1;
            }
        }
    }

    fn set_ref_alleles(&self, m: i32, input_alleles: &mut [i32], al_cnts: &mut [i32]) {
        let ref_gt = self.opt_ref.as_ref().expect("ref present");
        let mut ref_hap = 0;
        let mut h1 = self.n_targ_haps;
        while h1 < self.n_haps {
            let h2 = h1 | 0b1;
            let a1 = ref_gt.allele(m, ref_hap);
            ref_hap += 1;
            let a2 = ref_gt.allele(m, ref_hap);
            ref_hap += 1;
            input_alleles[h1 as usize] = a1;
            input_alleles[h2 as usize] = a2;
            al_cnts[a1 as usize] += 1;
            al_cnts[a2 as usize] += 1;
            h1 += 2;
        }
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
        let gt = std::env::temp_dir().join("beagle_rs_prp_gt.vcf");
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
    fn first_marker_copies_input_alleles_and_cdf() {
        let fpd = fpd();
        let mut p = PbwtRecPhaser::new(&fpd);
        let n_haps = p.n_haps();
        assert_eq!(n_haps, 6); // 3 diploid samples, no ref
        let mut alleles = vec![-1i32; n_haps as usize];
        let mut missing = vec![false; 3];
        let mut unph_het = vec![false; 3];
        let al_cnts = p.phase(-1, &mut alleles, 0, &mut missing, &mut unph_het);
        // marker0: S0=0|0, S1=0|1, S2=1|1 -> haps [0,0,0,1,1,1]
        assert_eq!(alleles, vec![0, 0, 0, 1, 1, 1]);
        // allele-count CDF: allele0 count 3, allele1 count 3 -> CDF [3,6]
        assert_eq!(al_cnts, vec![3, 6]);
        // S1 is an unphased het at the first marker (no overlap) but gets phased by PBWT
        assert!(!missing.iter().any(|&m| m));
    }
}
