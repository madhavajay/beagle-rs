//! Port of `phase/PbwtPhaser.java` — produces the initial per-sample phasing by running a
//! `FwdPbwtPhaser` over each high-frequency cM window and stitching the windows together
//! with phase alignment across their overlaps.

use crate::blbutil::DoubleArray;
use crate::ints::{IntArray, WrappedIntArray};

use super::{FixedPhaseData, FwdPbwtPhaser, SamplePhase};

/// Port of `phase/PbwtPhaser.java`.
pub struct PbwtPhaser {
    start: i32,
    end: i32,
    fwd_pbwt: FwdPbwtPhaser,
}

struct Indices {
    miss_indices: WrappedIntArray,
    het_indices: WrappedIntArray,
}

fn ins_pt(list: &WrappedIntArray, value: i32) -> i32 {
    let index = list.binary_search(value);
    if index < 0 {
        -index - 1
    } else {
        index
    }
}

fn from(gen_pos: &DoubleArray, pos: f64) -> i32 {
    let ins = gen_pos.binary_search(pos);
    if ins < 0 {
        -ins - 1
    } else {
        ins
    }
}

fn to(gen_pos: &DoubleArray, pos: f64) -> i32 {
    let ins = gen_pos.binary_search(pos);
    if ins < 0 {
        -ins - 1
    } else {
        ins + 1
    }
}

/// `alignmentHet` — the marker used to align this window's phase with the previous window,
/// or `-1` if none.
fn alignment_het(het_list: &WrappedIntArray, start: i32, copy_start: i32, overlap_end: i32) -> i32 {
    if het_list.size() == 0 {
        return -1;
    }
    let mut index = ins_pt(het_list, copy_start);
    if index == het_list.size() || (het_list.get(index) >= overlap_end && index > 0) {
        index -= 1;
    }
    let het = het_list.get(index);
    if start <= het && het < overlap_end {
        het
    } else {
        -1
    }
}

fn hi_freq_windows(fpd: &FixedPhaseData) -> Vec<[i32; 2]> {
    let gen_pos = fpd.stage1_map().gen_pos();
    let n_markers = gen_pos.size();
    let n_threads = fpd.par().nthreads();
    let total_cm = gen_pos.get(gen_pos.size() - 1) - gen_pos.get(0);
    let overlap_cm = 0.5f64;
    let advance_cm = (4.0 * overlap_cm).max(total_cm / n_threads as f64);
    let mut window_list = Vec::new();
    let mut frm = 0;
    let mut t = to(gen_pos, gen_pos.get(frm) + advance_cm);
    while t < n_markers {
        window_list.push([frm, t]);
        frm = from(gen_pos, gen_pos.get(t) - overlap_cm);
        t = to(gen_pos, gen_pos.get(t) + advance_cm);
    }
    debug_assert_eq!(t, n_markers);
    window_list.push([frm, t]);
    window_list
}

fn pbwt_phasers(fpd: &FixedPhaseData, seed: i64) -> Vec<PbwtPhaser> {
    let windows = hi_freq_windows(fpd);
    windows
        .iter()
        .enumerate()
        .map(|(j, w)| PbwtPhaser::new(fpd, w[0], w[1], seed + j as i64))
        .collect()
}

fn int_lists(length: i32) -> Vec<crate::ints::IntList> {
    (0..length).map(|_| crate::ints::IntList::new()).collect()
}

fn indices(fpd: &FixedPhaseData, s_start: i32, s_end: i32) -> Vec<Indices> {
    let gt = fpd.stage1_targ_gt();
    let overlap = fpd.stage1_overlap();
    let n_markers = gt.n_markers();
    let len = s_end - s_start;
    let mut miss_indices = int_lists(len);
    let mut het_indices = int_lists(len);
    let mut not_first_het = vec![false; len as usize];
    for m in 0..n_markers {
        for s in s_start..s_end {
            let ss = (s - s_start) as usize;
            let hap1 = s << 1;
            let a1 = gt.allele(m, hap1);
            let a2 = gt.allele(m, hap1 | 0b1);
            if a1 < 0 || a2 < 0 {
                miss_indices[ss].add(m);
            } else if a1 != a2 {
                if m >= overlap && not_first_het[ss] {
                    het_indices[ss].add(m);
                } else {
                    not_first_het[ss] = true;
                }
            }
        }
    }
    (0..len as usize)
        .map(|j| Indices {
            miss_indices: WrappedIntArray::from_list(&miss_indices[j]),
            het_indices: WrappedIntArray::from_list(&het_indices[j]),
        })
        .collect()
}

impl PbwtPhaser {
    fn new(fpd: &FixedPhaseData, start: i32, end: i32, seed: i64) -> Self {
        assert!(
            start >= 0 && end <= fpd.targ_gt().n_markers() && start < end,
            "{start}"
        );
        PbwtPhaser {
            start,
            end,
            fwd_pbwt: FwdPbwtPhaser::new(fpd, start, end, seed),
        }
    }

    /// `PbwtPhaser.initPhase(FixedPhaseData fpd, long seed)` — the initial per-sample phasing.
    pub fn init_phase(fpd: &FixedPhaseData, seed: i64) -> Vec<SamplePhase> {
        let ppa = pbwt_phasers(fpd, seed);
        let n_samples = fpd.stage1_targ_gt().n_samples();
        let n_threads = fpd.par().nthreads();
        let max_step_size = 128;
        let step_size = ((n_samples + n_threads - 1) / n_threads).min(max_step_size);
        let n_steps = (n_samples + (step_size - 1)) / step_size;
        let mut phase: Vec<Option<SamplePhase>> = (0..n_samples).map(|_| None).collect();
        for step in 0..n_steps {
            set_sample_phase(fpd, &ppa, &mut phase, step, step_size);
        }
        phase.into_iter().map(|p| p.expect("phase")).collect()
    }

    fn switch_hap_labels(&self, sample: i32, hap1: &[i32], hap2: &[i32], align_het: i32) -> bool {
        let h1 = sample << 1;
        let h2 = h1 | 0b1;
        let a1 = hap1[align_het as usize];
        let a2 = hap2[align_het as usize];
        let b1 = self.fwd_pbwt.allele(align_het, h1);
        let b2 = self.fwd_pbwt.allele(align_het, h2);
        a1 == b2 && a2 == b1
    }

    fn copy_haps(
        &self,
        haps: &mut [Vec<i32>],
        indices: &[Indices],
        overlap_end: i32,
        s_start: i32,
        s_end: i32,
    ) {
        let copy_start = (self.start + overlap_end) >> 1;
        let len = (s_end - s_start) as usize;
        let mut switched = vec![false; len];
        if self.start > 0 {
            for s in s_start..s_end {
                let ss = (s - s_start) as usize;
                let hh1 = ss << 1;
                let hh2 = hh1 | 0b1;
                let align_het = alignment_het(
                    &indices[ss].het_indices,
                    self.start,
                    copy_start,
                    overlap_end,
                );
                if align_het >= 0 && self.switch_hap_labels(s, &haps[hh1], &haps[hh2], align_het) {
                    switched[ss] = true;
                }
            }
        }
        for s in s_start..s_end {
            let ss = (s - s_start) as usize;
            let hh1 = ss << 1;
            let hh2 = hh1 | 0b1;
            let h1 = s << 1;
            let h2 = h1 | 0b1;
            let (d1, d2) = if switched[ss] { (hh2, hh1) } else { (hh1, hh2) };
            for m in copy_start..self.end {
                haps[d1][m as usize] = self.fwd_pbwt.allele(m, h1);
                haps[d2][m as usize] = self.fwd_pbwt.allele(m, h2);
            }
        }
    }
}

fn set_sample_phase(
    fpd: &FixedPhaseData,
    ppa: &[PbwtPhaser],
    phase: &mut [Option<SamplePhase>],
    step: i32,
    step_size: i32,
) {
    let gt = fpd.stage1_targ_gt();
    let s_start = step * step_size;
    let s_end = (s_start + step_size).min(gt.n_samples());
    let idx = indices(fpd, s_start, s_end);

    let mut haps: Vec<Vec<i32>> =
        vec![vec![0i32; gt.n_markers() as usize]; ((s_end - s_start) << 1) as usize];
    ppa[0].copy_haps(&mut haps, &idx, 0, s_start, s_end);
    for j in 1..ppa.len() {
        let overlap_end = ppa[j - 1].end;
        ppa[j].copy_haps(&mut haps, &idx, overlap_end, s_start, s_end);
    }

    for s in s_start..s_end {
        let ss = (s - s_start) as usize;
        let hh1 = ss << 1;
        let hh2 = hh1 | 0b1;
        phase[s as usize] = Some(SamplePhase::new(
            s,
            gt.markers().clone(),
            fpd.stage1_map().gen_pos(),
            &haps[hh1],
            &haps[hh2],
            &idx[ss].het_indices,
            &idx[ss].miss_indices,
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, GT, HEADER_PREFIX,
    };
    use std::io::Write;
    use std::rc::Rc;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_pp_gt.vcf");
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
    fn init_phase_produces_one_phasing_per_sample() {
        let fpd = fpd();
        let phase = PbwtPhaser::init_phase(&fpd, 99999);
        assert_eq!(phase.len(), 3);
        let n = fpd.stage1_targ_gt().n_markers();
        for m in 0..n {
            // S0 = 0|0
            assert_eq!(phase[0].allele1(m), 0);
            assert_eq!(phase[0].allele2(m), 0);
            // S1 het (some phase)
            let (a, b) = (phase[1].allele1(m), phase[1].allele2(m));
            assert!((a, b) == (0, 1) || (a, b) == (1, 0));
            // S2 = 1|1
            assert_eq!(phase[2].allele1(m), 1);
            assert_eq!(phase[2].allele2(m), 1);
        }
    }
}
