//! Port of `phase/LowFreqPbwtPhaseIbs.java` — like `PbwtPhaseIbs`, but augments the long-IBS
//! PBWT search with rare-variant ("stage-2") sharing: low-frequency alleles in a step pin
//! candidate haplotypes together, and the nearest such candidate within a backoff window
//! (skipping IBS2 matches) is preferred over the random long-IBS match.

use crate::beagleutil::PbwtDivUpdater;
use crate::ints::{IntArray, IntList, WrappedIntArray};
use crate::jdk::Random;
use crate::vcf::{Steps, XRefGT};

use super::{CodedSteps, FixedPhaseData, Ibs2, PbwtIbsData, PhaseData};

/// Port of `phase/LowFreqPbwtPhaseIbs.java`.
pub struct LowFreqPbwtPhaseIbs<'a> {
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

impl<'a> LowFreqPbwtPhaseIbs<'a> {
    /// `new LowFreqPbwtPhaseIbs(PhaseData, CodedSteps, boolean useBwd)`.
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
        LowFreqPbwtPhaseIbs {
            phase_data,
            all_haps,
            ibs_haps,
        }
    }

    /// `phaseData()`.
    pub fn phase_data(&self) -> &PhaseData {
        self.phase_data
    }

    /// `allHaps()`.
    pub fn all_haps(&self) -> &XRefGT {
        &self.all_haps
    }

    /// `ibsHap(int hap, int step)`.
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
    let mut a_inv = vec![0i32; n_haps as usize];
    let mut i_to_prev_i = vec![0i32; n_haps as usize];
    let mut i_to_next_i = vec![0i32; n_haps as usize];
    let mut ibs_haps0: Vec<Option<WrappedIntArray>> =
        (0..(end_step - start_step)).map(|_| None).collect();

    for j in (end_step..buffer_end_step).rev() {
        let ia = data.coded_steps().get(j);
        pbwt.bwd_update(ia.int_array(), ia.value_size(), j, &mut a, &mut d);
    }
    for j in (start_step..end_step).rev() {
        let ia = data.coded_steps().get(j);
        pbwt.bwd_update(ia.int_array(), ia.value_size(), j, &mut a, &mut d);
        set_inv(&a, &mut a_inv);
        set_i_to_prev_next_i(
            phase_data,
            data.coded_steps().steps(),
            j,
            &a_inv,
            &mut i_to_prev_i,
            &mut i_to_next_i,
        );
        ibs_haps0[(j - start_step) as usize] = Some(get_bwd_ibs_haps(
            phase_data,
            j,
            &a,
            &mut d,
            &i_to_prev_i,
            &i_to_next_i,
            data,
        ));
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
    let mut a_inv = vec![0i32; n_haps as usize];
    let mut i_to_prev_i = vec![0i32; n_haps as usize];
    let mut i_to_next_i = vec![0i32; n_haps as usize];
    let mut ibs_haps0: Vec<Option<WrappedIntArray>> =
        (0..(end_step - start_step)).map(|_| None).collect();

    for j in buffer_start_step..start_step {
        let ia = data.coded_steps().get(j);
        pbwt.fwd_update(ia.int_array(), ia.value_size(), j, &mut a, &mut d);
    }
    for j in start_step..end_step {
        let ia = data.coded_steps().get(j);
        pbwt.fwd_update(ia.int_array(), ia.value_size(), j, &mut a, &mut d);
        set_inv(&a, &mut a_inv);
        set_i_to_prev_next_i(
            phase_data,
            data.coded_steps().steps(),
            j,
            &a_inv,
            &mut i_to_prev_i,
            &mut i_to_next_i,
        );
        ibs_haps0[(j - start_step) as usize] = Some(get_fwd_ibs_haps(
            phase_data,
            j,
            &a,
            &mut d,
            &i_to_prev_i,
            &i_to_next_i,
            data,
        ));
    }
    ibs_haps0.into_iter().map(|h| h.expect("ibs")).collect()
}

#[allow(clippy::too_many_arguments)]
fn get_bwd_ibs_haps(
    phase_data: &PhaseData,
    step: i32,
    a: &[i32],
    d: &mut [i32],
    i_to_prev_i: &[i32],
    i_to_next_i: &[i32],
    data: &PbwtIbsData,
) -> WrappedIntArray {
    let mut rand = Random::new(phase_data.seed() + step as i64);
    let ibs2 = phase_data.fpd().stage1_ibs2();
    let m_start = data.coded_steps().steps().start(step);
    let m_incl_end = data.coded_steps().steps().end(step) - 1;
    let n_targ_haps = data.n_targ_haps();
    let n_candidates = data.n_candidates();
    let mut selected_haps = vec![0i32; n_targ_haps as usize];
    let n = a.len();
    d[n] = step - 1; // set sentinel
    for i in 0..n {
        if a[i] < n_targ_haps {
            let best_i = best_bwd_stage2_index(
                ibs2,
                step,
                m_start,
                m_incl_end,
                i as i32,
                a,
                d,
                i_to_prev_i,
                i_to_next_i,
                data,
            );
            if best_i >= 0 {
                selected_haps[a[i] as usize] = a[best_i as usize];
            } else {
                let mut u = i as i32;
                let mut v = i as i32 + 1;
                let mut u_next_match_end = d[u as usize];
                let mut v_next_match_end = d[v as usize];
                while (v - u) < n_candidates
                    && (step <= u_next_match_end || step <= v_next_match_end)
                {
                    if u_next_match_end <= v_next_match_end {
                        v += 1;
                        v_next_match_end = d[v as usize].min(v_next_match_end);
                    } else {
                        u -= 1;
                        u_next_match_end = d[u as usize].min(u_next_match_end);
                    }
                }
                selected_haps[a[i] as usize] =
                    get_match(ibs2, m_start, m_incl_end, i as i32, u, v, a, &mut rand);
            }
        }
    }
    WrappedIntArray::from_slice(&selected_haps)
}

#[allow(clippy::too_many_arguments)]
fn get_fwd_ibs_haps(
    phase_data: &PhaseData,
    step: i32,
    a: &[i32],
    d: &mut [i32],
    i_to_prev_i: &[i32],
    i_to_next_i: &[i32],
    data: &PbwtIbsData,
) -> WrappedIntArray {
    let mut rand = Random::new(phase_data.seed() + step as i64);
    let ibs2 = phase_data.fpd().stage1_ibs2();
    let m_start = data.coded_steps().steps().start(step);
    let m_incl_end = data.coded_steps().steps().end(step) - 1;
    let n_targ_haps = data.n_targ_haps();
    let n_candidates = data.n_candidates();
    let mut selected_haps = vec![0i32; n_targ_haps as usize];
    let n = a.len();
    d[n] = step + 1; // set sentinel
    for i in 0..n {
        if a[i] < n_targ_haps {
            let best_i = best_fwd_stage2_index(
                ibs2,
                step,
                m_start,
                m_incl_end,
                i as i32,
                a,
                d,
                i_to_prev_i,
                i_to_next_i,
                data,
            );
            if best_i >= 0 {
                selected_haps[a[i] as usize] = a[best_i as usize];
            } else {
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
                selected_haps[a[i] as usize] =
                    get_match(ibs2, m_start, m_incl_end, i as i32, u, v, a, &mut rand);
            }
        }
    }
    WrappedIntArray::from_slice(&selected_haps)
}

#[allow(clippy::too_many_arguments)]
fn best_fwd_stage2_index(
    ibs2: &Ibs2,
    step: i32,
    m_start: i32,
    m_incl_end: i32,
    i: i32,
    a: &[i32],
    d: &[i32],
    i_to_prev_i: &[i32],
    i_to_next_i: &[i32],
    tmp_data: &PbwtIbsData,
) -> i32 {
    let n = a.len() as i32;
    let mut best_prev_match = -1;
    let mut best_next_match = -1;
    let mut prev_match_start = 0;
    let mut next_match_start = 0;

    let min_match_start = if (i + 1) < n {
        d[i as usize].min(d[(i + 1) as usize])
    } else {
        d[i as usize]
    };
    let d_max = (min_match_start + tmp_data.max_backoff_steps()).min(step);
    let mut prev_i = i_to_prev_i[i as usize];
    while prev_i > i32::MIN
        && ibs2.are_ibs2_in(
            a[i as usize] >> 1,
            a[prev_i as usize] >> 1,
            m_start,
            m_incl_end,
        )
    {
        prev_i = i_to_prev_i[prev_i as usize];
    }
    if prev_i > i32::MIN {
        debug_assert!(prev_i < i);
        let mut u = i;
        while (u - 1) != prev_i && d[u as usize] <= d_max {
            prev_match_start = prev_match_start.max(d[u as usize]);
            u -= 1;
        }
        if (u - 1) == prev_i && d[u as usize] <= d_max {
            prev_match_start = prev_match_start.max(d[u as usize]);
            best_prev_match = prev_i;
        }
    }
    let mut next_i = i_to_next_i[i as usize];
    while next_i < i32::MAX
        && ibs2.are_ibs2_in(
            a[i as usize] >> 1,
            a[next_i as usize] >> 1,
            m_start,
            m_incl_end,
        )
    {
        next_i = i_to_next_i[next_i as usize];
    }
    if next_i < i32::MAX {
        debug_assert!(i < next_i);
        let mut v = i;
        while (v + 1) != next_i && d[(v + 1) as usize] <= d_max {
            v += 1;
            next_match_start = next_match_start.max(d[v as usize]);
        }
        if (v + 1) == next_i && d[(v + 1) as usize] <= d_max {
            v += 1;
            next_match_start = next_match_start.max(d[v as usize]);
            best_next_match = next_i;
        }
    }
    if prev_match_start < next_match_start && best_prev_match != -1 {
        best_prev_match
    } else {
        best_next_match
    }
}

#[allow(clippy::too_many_arguments)]
fn best_bwd_stage2_index(
    ibs2: &Ibs2,
    step: i32,
    m_start: i32,
    m_incl_end: i32,
    i: i32,
    a: &[i32],
    d: &[i32],
    i_to_prev_i: &[i32],
    i_to_next_i: &[i32],
    data: &PbwtIbsData,
) -> i32 {
    let n = a.len() as i32;
    let n_steps_m1 = data.coded_steps().steps().size() - 1;
    let mut best_prev_match = -1;
    let mut best_next_match = -1;
    let mut prev_match_incl_end = n_steps_m1;
    let mut next_match_incl_end = n_steps_m1;
    let max_match_start = if (i + 1) < n {
        d[i as usize].max(d[(i + 1) as usize])
    } else {
        d[i as usize]
    };
    let d_min = (max_match_start - data.max_backoff_steps()).max(step);
    let mut prev_i = i_to_prev_i[i as usize];
    while prev_i > i32::MIN
        && ibs2.are_ibs2_in(
            a[i as usize] >> 1,
            a[prev_i as usize] >> 1,
            m_start,
            m_incl_end,
        )
    {
        prev_i = i_to_prev_i[prev_i as usize];
    }
    if prev_i > i32::MIN {
        debug_assert!(prev_i < i);
        let mut u = i;
        while (u - 1) != prev_i && d[u as usize] >= d_min {
            prev_match_incl_end = prev_match_incl_end.min(d[u as usize]);
            u -= 1;
        }
        if (u - 1) == prev_i && d[u as usize] >= d_min {
            prev_match_incl_end = prev_match_incl_end.min(d[u as usize]);
            best_prev_match = prev_i;
        }
    }
    let mut next_i = i_to_next_i[i as usize];
    while next_i < i32::MAX
        && ibs2.are_ibs2_in(
            a[i as usize] >> 1,
            a[next_i as usize] >> 1,
            m_start,
            m_incl_end,
        )
    {
        next_i = i_to_next_i[next_i as usize];
    }
    if next_i < i32::MAX {
        debug_assert!(i < next_i);
        let mut v = i;
        while (v + 1) != next_i && d[(v + 1) as usize] >= d_min {
            v += 1;
            next_match_incl_end = next_match_incl_end.min(d[v as usize]);
        }
        if (v + 1) == next_i && d[(v + 1) as usize] >= d_min {
            v += 1;
            next_match_incl_end = next_match_incl_end.min(d[v as usize]);
            best_next_match = next_i;
        }
    }
    if prev_match_incl_end > next_match_incl_end && best_prev_match != -1 {
        best_prev_match
    } else {
        best_next_match
    }
}

#[allow(clippy::too_many_arguments)]
fn get_match(
    ibs2: &Ibs2,
    m_start: i32,
    m_incl_end: i32,
    i: i32,
    i_start: i32,
    i_end: i32,
    a: &[i32],
    rand: &mut Random,
) -> i32 {
    let i_length = i_end - i_start;
    if i_length == 1 {
        return -1;
    }
    let mut found = -1;
    let mut index = i_start + rand.next_int_bound(i_length);
    let mut j = 0;
    while j < i_length && found == -1 {
        if !ibs2.are_ibs2_in(
            a[i as usize] >> 1,
            a[index as usize] >> 1,
            m_start,
            m_incl_end,
        ) {
            found = a[index as usize];
        }
        index += 1;
        if index == i_end {
            index = i_start;
        }
        j += 1;
    }
    found
}

fn set_i_to_prev_next_i(
    phase_data: &PhaseData,
    steps: &Steps,
    step: i32,
    inv_a: &[i32],
    i_to_prev_i: &mut [i32],
    i_to_next_i: &mut [i32],
) {
    i_to_prev_i.iter_mut().for_each(|x| *x = i32::MIN);
    i_to_next_i.iter_mut().for_each(|x| *x = i32::MAX);
    let low_freq_hap_lists = low_freq_hap_lists(phase_data, steps, step);
    for haps in &low_freq_hap_lists {
        let idx = sorted_a_indices(haps, inv_a);
        for k in 1..idx.len() {
            let i0 = idx[k - 1];
            let i1 = idx[k];
            if i0 > i_to_prev_i[i1 as usize] {
                i_to_prev_i[i1 as usize] = i0;
            }
            if i1 < i_to_next_i[i0 as usize] {
                i_to_next_i[i0 as usize] = i1;
            }
        }
    }
}

fn sorted_a_indices(haps: &IntList, inv_a: &[i32]) -> Vec<i32> {
    let mut out: Vec<i32> = (0..haps.size())
        .map(|j| inv_a[haps.get(j) as usize])
        .collect();
    out.sort_unstable();
    out
}

fn low_freq_hap_lists(phase_data: &PhaseData, steps: &Steps, step: i32) -> Vec<IntList> {
    let fpd = phase_data.fpd();
    let hi_freq_indices = fpd.stage1_to2();
    let start = if step == 0 {
        0
    } else {
        hi_freq_indices.get(steps.start(step))
    };
    let end = if step + 1 < steps.size() {
        hi_freq_indices.get(steps.start(step + 1))
    } else {
        fpd.targ_gt().n_markers()
    };
    low_freq_hap_lists_range(fpd, start, end)
}

fn low_freq_hap_lists_range(fpd: &FixedPhaseData, start: i32, end: i32) -> Vec<IntList> {
    let mut hap_lists = Vec::new();
    let markers = fpd.targ_gt().markers().clone();
    for m in start..end {
        let n_alleles = markers.marker(m).n_alleles();
        for al in 0..n_alleles {
            let carriers = fpd.carriers(m, al);
            if carriers.size() > 1 {
                hap_lists.push(hap_list(carriers.size(), |j| carriers.get(j)));
            }
        }
    }
    hap_lists
}

fn hap_list(n_carriers: i32, carrier: impl Fn(i32) -> i32) -> IntList {
    let mut hap_list = IntList::with_capacity(2 * n_carriers);
    for j in 0..n_carriers {
        let sample = carrier(j);
        let h1 = sample << 1;
        hap_list.add(h1);
        hap_list.add(h1 | 0b1);
    }
    hap_list
}

fn set_inv(a: &[i32], a_inv: &mut [i32]) {
    for (j, &aj) in a.iter().enumerate() {
        a_inv[aj as usize] = j as i32;
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
        let gt = std::env::temp_dir().join("beagle_rs_lfppi_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..6 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 6],
        );
        let mp = MarkerParser::new(true, true, true, true);
        // 8 markers with a mix; rare alt alleles in some markers exercise stage-2 sharing.
        let gts = [
            "0|0\t0|1\t1|1\t0|1\t0|0\t1|0",
            "0|1\t0|0\t1|0\t1|1\t0|1\t0|0",
            "1|1\t0|1\t0|0\t0|1\t1|0\t0|1",
            "0|0\t1|1\t0|1\t1|0\t0|0\t1|1",
            "0|1\t1|0\t1|1\t0|0\t0|1\t1|0",
            "1|0\t0|1\t0|0\t1|1\t1|0\t0|1",
            "0|0\t0|0\t0|1\t0|0\t0|0\t0|0", // rare alt at one hap
            "1|1\t1|1\t1|0\t1|1\t1|1\t1|1", // rare ref at one hap
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
            let ibs = LowFreqPbwtPhaseIbs::new(&pd, pd.coded_steps(), use_bwd);
            assert_eq!(ibs.all_haps().n_haps(), n_haps);
            for step in 0..n_steps {
                for hap in 0..n_haps {
                    let m = ibs.ibs_hap(hap, step);
                    assert!(m == -1 || (0..n_haps).contains(&m));
                    if m >= 0 {
                        assert_ne!(m, hap);
                    }
                }
            }
        }
    }
}
