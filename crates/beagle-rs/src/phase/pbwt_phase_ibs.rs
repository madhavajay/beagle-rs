//! Port of `phase/PbwtPhaseIbs.java` — uses the PBWT (Durbin 2014) to find, for each target
//! haplotype and each step, a long IBS haplotype that spans the step (excluding IBS2 matches).
//!
//! Java runs the batches as a parallel stream and flat-maps the per-batch results in batch
//! order; the batches partition the steps contiguously, so the sequential build produces the
//! identical `ibsHaps[step]` array.

use crate::beagleutil::PbwtDivUpdater;
use crate::ints::{IntArray, WrappedIntArray};
use crate::jdk::Random;
use crate::vcf::XRefGT;

use super::{CodedSteps, PbwtIbsData, PhaseData};

/// Port of `phase/PbwtPhaseIbs.java`.
pub struct PbwtPhaseIbs<'a> {
    phase_data: &'a PhaseData,
    all_haps: XRefGT,
    ibs_haps: Vec<WrappedIntArray>, // [step] -> selected IBS hap per target hap
}

fn check_consistency(phase_data: &PhaseData, coded_steps: &CodedSteps) {
    let fpd = phase_data.fpd();
    let consistent = fpd.stage1_steps().size() == coded_steps.steps().size()
        && fpd.stage1_xref_gt().is_some() == coded_steps.ref_haps().is_some()
        && fpd.targ_gt().samples() == coded_steps.targ_samples();
    assert!(consistent, "inconsistent data");
}

impl<'a> PbwtPhaseIbs<'a> {
    /// `new PbwtPhaseIbs(PhaseData, CodedSteps, boolean useBwd)`.
    pub fn new(phase_data: &'a PhaseData, coded_steps: CodedSteps, use_bwd: bool) -> Self {
        check_consistency(phase_data, &coded_steps);
        let all_haps = coded_steps.all_haps().clone();
        let data = PbwtIbsData::new(phase_data, coded_steps);
        let mut ibs_haps: Vec<WrappedIntArray> = Vec::with_capacity(data.n_steps() as usize);
        for batch in 0..data.n_batches() {
            let batch_haps = if use_bwd {
                bwd_ibs_haps(phase_data, &data, batch)
            } else {
                fwd_ibs_haps(phase_data, &data, batch)
            };
            ibs_haps.extend(batch_haps);
        }
        PbwtPhaseIbs {
            phase_data,
            all_haps,
            ibs_haps,
        }
    }

    /// `phaseData()`.
    pub fn phase_data(&self) -> &PhaseData {
        self.phase_data
    }

    /// `allHaps()` — phased target + reference genotypes.
    pub fn all_haps(&self) -> &XRefGT {
        &self.all_haps
    }

    /// `ibsHap(int hap, int step)` — an IBS haplotype index for `hap` over `step`, or `-1`.
    pub fn ibs_hap(&self, hap: i32, step: i32) -> i32 {
        self.ibs_haps[step as usize].get(hap)
    }
}

fn bwd_ibs_haps(phase_data: &PhaseData, data: &PbwtIbsData, batch: i32) -> Vec<WrappedIntArray> {
    let start_step = data.start_step(batch);
    let end_step = data.end_step(batch);
    let buffer_end_step = data.buffer_end_step(end_step);

    let n_haps = data.n_haps();
    let mut pbwt = PbwtDivUpdater::new(n_haps);
    let mut a: Vec<i32> = (0..n_haps).collect();
    let mut d: Vec<i32> = vec![buffer_end_step - 1; (n_haps + 1) as usize]; // last entry is a sentinel
    let mut ibs_haps0: Vec<Option<WrappedIntArray>> =
        (0..(end_step - start_step)).map(|_| None).collect();

    for j in (end_step..buffer_end_step).rev() {
        let ia = data.coded_steps().get(j);
        pbwt.bwd_update(ia.int_array(), ia.value_size(), j, &mut a, &mut d);
    }
    for j in (start_step..end_step).rev() {
        let ia = data.coded_steps().get(j);
        pbwt.bwd_update(ia.int_array(), ia.value_size(), j, &mut a, &mut d);
        ibs_haps0[(j - start_step) as usize] =
            Some(get_bwd_ibs_haps(phase_data, j, &a, &mut d, data));
    }
    ibs_haps0.into_iter().map(|h| h.expect("ibs")).collect()
}

fn fwd_ibs_haps(phase_data: &PhaseData, data: &PbwtIbsData, batch: i32) -> Vec<WrappedIntArray> {
    let start_step = data.start_step(batch);
    let end_step = data.end_step(batch);
    let buffer_start_step = data.buffer_start_step(start_step);

    let n_haps = data.n_haps();
    let mut pbwt = PbwtDivUpdater::new(n_haps);
    let mut a: Vec<i32> = (0..n_haps).collect();
    let mut d: Vec<i32> = vec![buffer_start_step; (n_haps + 1) as usize]; // last entry is a sentinel
    let mut ibs_haps0: Vec<Option<WrappedIntArray>> =
        (0..(end_step - start_step)).map(|_| None).collect();

    for j in buffer_start_step..start_step {
        let ia = data.coded_steps().get(j);
        pbwt.fwd_update(ia.int_array(), ia.value_size(), j, &mut a, &mut d);
    }
    for j in start_step..end_step {
        let ia = data.coded_steps().get(j);
        pbwt.fwd_update(ia.int_array(), ia.value_size(), j, &mut a, &mut d);
        ibs_haps0[(j - start_step) as usize] =
            Some(get_fwd_ibs_haps(phase_data, j, &a, &mut d, data));
    }
    ibs_haps0.into_iter().map(|h| h.expect("ibs")).collect()
}

fn get_bwd_ibs_haps(
    phase_data: &PhaseData,
    step: i32,
    a: &[i32],
    d: &mut [i32],
    data: &PbwtIbsData,
) -> WrappedIntArray {
    let mut rand = Random::new(phase_data.seed() + step as i64);
    let m_start = data.coded_steps().steps().start(step);
    let m_incl_end = data.coded_steps().steps().end(step) - 1;
    let n_targ_haps = data.n_targ_haps();
    let n_candidates = data.n_candidates();
    let mut selected_haps = vec![0i32; n_targ_haps as usize];
    let ibs2 = phase_data.fpd().stage1_ibs2();
    let n = a.len();
    d[0] = step - 2; // set sentinels (old d[0], d[a.len()] need not be restored)
    d[n] = step - 2;
    for i in 0..n {
        if a[i] < n_targ_haps {
            let hap = a[i];
            let s1 = hap >> 1;
            let mut u = i as i32;
            let mut v = i as i32 + 1;
            let mut u_next_match_end = d[u as usize];
            let mut v_next_match_end = d[v as usize];
            while (v - u) < n_candidates && (step <= u_next_match_end || step <= v_next_match_end) {
                if u_next_match_end <= v_next_match_end {
                    v += 1;
                    v_next_match_end = d[v as usize].min(v_next_match_end);
                } else {
                    u -= 1;
                    u_next_match_end = d[u as usize].min(u_next_match_end);
                }
            }
            let count = v - u;
            selected_haps[hap as usize] = -1;
            if count > 1 {
                let mut index = u + rand.next_int_bound(count);
                for _ in 0..count {
                    if index == v {
                        index = u;
                    }
                    if index != i as i32
                        && !ibs2.are_ibs2_in(s1, a[index as usize] >> 1, m_start, m_incl_end)
                    {
                        selected_haps[hap as usize] = a[index as usize];
                        break;
                    }
                    index += 1;
                }
            }
        }
    }
    WrappedIntArray::from_slice(&selected_haps)
}

fn get_fwd_ibs_haps(
    phase_data: &PhaseData,
    step: i32,
    a: &[i32],
    d: &mut [i32],
    data: &PbwtIbsData,
) -> WrappedIntArray {
    let steps = phase_data.fpd().stage1_steps();
    let mut rand = Random::new(phase_data.seed() + step as i64);
    let n_targ_haps = phase_data.fpd().targ_gt().n_haps();
    let n_candidates = data.n_candidates();
    let m_start = steps.start(step);
    let m_incl_end = steps.end(step) - 1;
    let mut selected_haps = vec![0i32; n_targ_haps as usize];
    let ibs2 = phase_data.fpd().stage1_ibs2();
    let n = a.len();
    d[0] = step + 2; // set sentinels
    d[n] = step + 2;
    for i in 0..n {
        if a[i] < n_targ_haps {
            let hap = a[i];
            let s1 = hap >> 1;
            let mut u = i as i32;
            let mut v = i as i32 + 1;
            let mut u_next_match_start = d[u as usize];
            let mut v_next_match_start = d[v as usize];
            while (v - u) < n_candidates
                && (u_next_match_start <= step || v_next_match_start <= step)
            {
                if v_next_match_start <= u_next_match_start {
                    v += 1;
                    v_next_match_start = d[v as usize].max(v_next_match_start);
                } else {
                    u -= 1;
                    u_next_match_start = d[u as usize].max(u_next_match_start);
                }
            }
            let count = v - u;
            selected_haps[hap as usize] = -1;
            if count > 1 {
                let mut index = u + rand.next_int_bound(count);
                for _ in 0..count {
                    if index == v {
                        index = u;
                    }
                    if index != i as i32
                        && !ibs2.are_ibs2_in(s1, a[index as usize] >> 1, m_start, m_incl_end)
                    {
                        selected_haps[hap as usize] = a[index as usize];
                        break;
                    }
                    index += 1;
                }
            }
        }
    }
    WrappedIntArray::from_slice(&selected_haps)
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

    use super::super::FixedPhaseData;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_ppi_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..4 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 4],
        );
        let mp = MarkerParser::new(true, true, true, true);
        // 6 markers; mix of homs and hets so distinct haplotype sequences exist.
        let gts = [
            "0|0\t0|1\t1|1\t0|1",
            "0|1\t0|0\t1|0\t1|1",
            "1|1\t0|1\t0|0\t0|1",
            "0|0\t1|1\t0|1\t1|0",
            "0|1\t1|0\t1|1\t0|0",
            "1|0\t0|1\t0|0\t1|1",
        ];
        let recs: Vec<Rc<dyn GTRec>> = gts
            .iter()
            .enumerate()
            .map(|(i, g)| {
                let pos = 1_000_000 + i as i32 * 1_000_000;
                let line = format!("chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t{g}");
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
    fn fwd_and_bwd_select_valid_ibs_haps() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let n_haps = pd.fpd().targ_gt().n_haps();
        let n_steps = pd.fpd().stage1_steps().size();

        for &use_bwd in &[false, true] {
            let ibs = PbwtPhaseIbs::new(&pd, pd.coded_steps(), use_bwd);
            assert_eq!(ibs.all_haps().n_haps(), n_haps);
            for step in 0..n_steps {
                for hap in 0..n_haps {
                    let m = ibs.ibs_hap(hap, step);
                    // selected hap is either -1 (none) or a valid haplotype != the query hap
                    assert!(m == -1 || (0..n_haps).contains(&m));
                    if m >= 0 {
                        assert_ne!(m, hap);
                    }
                }
            }
        }
    }
}
