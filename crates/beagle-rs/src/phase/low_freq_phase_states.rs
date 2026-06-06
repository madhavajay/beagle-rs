//! Port of `phase/LowFreqPhaseStates.java` — builds the Li & Stephens HMM states for a single
//! target *haplotype*, enriched for reference haplotypes carrying low-frequency variants. Each
//! composite reference haplotype is stored as a list of (hap, end-marker) segments and
//! evaluated lazily while copying the per-marker state alleles and mismatches.

use crate::beagleutil::CompHapSegment;
use crate::blbutil::Utilities;
use crate::ints::{IntIntMap, IntList};
use crate::jdk::{PriorityQueue, Random};
use crate::vcf::GT;

use super::LowFreqPhaseIbs;

const NIL: i32 = -103;

/// Port of `phase/LowFreqPhaseStates.java`.
pub struct LowFreqPhaseStates<'a> {
    ibs_haps: &'a LowFreqPhaseIbs<'a>,
    n_markers: i32,
    max_states: i32,
    min_steps: i32,
    hap_to_last_ibs_step: IntIntMap,
    q: PriorityQueue<CompHapSegment>,
    comp_hap_hap: Vec<IntList>,
    comp_hap_end: Vec<IntList>,
    segment_index: Vec<i32>,
    comp_hap_to_hap: Vec<i32>,
    comp_hap_to_end: Vec<i32>,
}

impl<'a> LowFreqPhaseStates<'a> {
    /// `new LowFreqPhaseStates(LowFreqPhaseIbs ibsHaps, int maxStates)`.
    pub fn new(ibs_haps: &'a LowFreqPhaseIbs<'a>, max_states: i32) -> Self {
        assert!(max_states >= 1, "{max_states}");
        let phase_data = ibs_haps.phase_data();
        let n_markers = ibs_haps.all_haps().n_markers();
        let phase_step = phase_data.fpd().ibs_step();
        let min_steps = 200.max((1.0f32 / phase_step).ceil() as i32);
        LowFreqPhaseStates {
            ibs_haps,
            n_markers,
            max_states,
            min_steps,
            hap_to_last_ibs_step: IntIntMap::new(max_states),
            q: PriorityQueue::new(max_states as usize),
            comp_hap_hap: (0..max_states).map(|_| IntList::new()).collect(),
            comp_hap_end: (0..max_states).map(|_| IntList::new()).collect(),
            segment_index: vec![0; max_states as usize],
            comp_hap_to_hap: vec![0; max_states as usize],
            comp_hap_to_end: vec![0; max_states as usize],
        }
    }

    /// `nTargHaps()`.
    pub fn n_targ_haps(&self) -> i32 {
        self.ibs_haps.phase_data().fpd().targ_gt().n_haps()
    }

    /// `nMarkers()`.
    pub fn n_markers(&self) -> i32 {
        self.ibs_haps.phase_data().fpd().targ_gt().n_markers()
    }

    /// `maxStates()`.
    pub fn max_states(&self) -> i32 {
        self.max_states
    }

    /// `ibsStates(int targHap, int[][] haps, byte[][] nMismatches)` — `haps[m][j]` is the
    /// `j`-th state's reference haplotype at marker `m`; `nMismatches[m][j]` is 0/1.
    pub fn ibs_states(
        &mut self,
        targ_hap: i32,
        haps: &mut [Vec<i32>],
        n_mismatches: &mut [Vec<u8>],
    ) -> i32 {
        let n_comp_haps = self.set_comp_ref_haps(targ_hap);
        self.copy_data(targ_hap, n_comp_haps, haps, n_mismatches);
        n_comp_haps
    }

    fn set_comp_ref_haps(&mut self, targ_hap: i32) -> i32 {
        let ibs_haps = self.ibs_haps;
        let steps = ibs_haps.phase_data().fpd().stage1_steps();
        self.q.clear();
        self.hap_to_last_ibs_step.clear();
        for j in 0..self.max_states as usize {
            self.comp_hap_hap[j].clear();
            self.comp_hap_end[j].clear();
        }
        let n = steps.size();
        for step in 0..n {
            self.add_ibs_hap(ibs_haps.fwd_ibs_hap(targ_hap, step), step);
            self.add_ibs_hap(ibs_haps.bwd_ibs_hap(targ_hap, step), step);
        }
        if self.q.is_empty() {
            self.fill_q_with_random_haps(targ_hap);
        }
        self.set_final_ref_segs()
    }

    fn add_ibs_hap(&mut self, ibs_hap: i32, step: i32) {
        if ibs_hap < 0 {
            return;
        }
        if self.hap_to_last_ibs_step.get(ibs_hap, NIL) == NIL {
            // hap is not currently in q
            self.update_head_of_q();
            let backoff = !self.q.is_empty()
                && (step - self.q.peek().expect("nonempty").last_ibs_step()) >= self.min_steps;
            if self.q.size() == self.max_states || backoff {
                let mut head = self.q.poll().expect("nonempty");
                let index = head.comp_hap_index() as usize;
                let prev_hap = head.hap();
                let next_start = self
                    .ibs_haps
                    .phase_data()
                    .fpd()
                    .stage1_steps()
                    .start((head.last_ibs_step() + step) >> 1);
                self.hap_to_last_ibs_step.remove(prev_hap);
                self.comp_hap_hap[index].add(ibs_hap); // hap of new segment
                self.comp_hap_end[index].add(next_start); // end of old segment
                head.update_segment(ibs_hap, next_start, step);
                self.q.add(head);
            } else {
                let index = self.q.size() as usize;
                self.comp_hap_hap[index].add(ibs_hap); // hap of new segment
                self.q
                    .add(CompHapSegment::new(ibs_hap, 0, step, index as i32));
            }
        }
        self.hap_to_last_ibs_step.put(ibs_hap, step);
    }

    fn update_head_of_q(&mut self) {
        if self.q.is_empty() {
            return;
        }
        loop {
            let (head_hap, head_last) = {
                let head = self.q.peek().expect("nonempty");
                (head.hap(), head.last_ibs_step())
            };
            let last_ibs_step = self.hap_to_last_ibs_step.get(head_hap, NIL);
            if head_last == last_ibs_step {
                break;
            }
            let mut head = self.q.poll().expect("nonempty");
            head.set_last_ibs_step(last_ibs_step);
            self.q.offer(head);
        }
    }

    fn set_final_ref_segs(&mut self) -> i32 {
        let n_comp_haps = self.q.size();
        let n_markers = self.n_markers;
        while let Some(head) = self.q.poll() {
            let comp_hap = head.comp_hap_index() as usize;
            self.comp_hap_end[comp_hap].add(n_markers); // add missing end of last segment
            self.segment_index[comp_hap] = 0;
            self.comp_hap_to_hap[comp_hap] = self.comp_hap_hap[comp_hap].get(0);
            self.comp_hap_to_end[comp_hap] = self.comp_hap_end[comp_hap].get(0);
        }
        n_comp_haps
    }

    // Indices `m`/`j` address parallel [marker][state] arrays plus the per-state segment
    // cursors; iterator form would obscure the lazy segment advance.
    #[allow(clippy::needless_range_loop)]
    fn copy_data(
        &mut self,
        targ_hap: i32,
        n_comp_haps: i32,
        haps: &mut [Vec<i32>],
        n_mismatches: &mut [Vec<u8>],
    ) {
        let all_haps = self.ibs_haps.all_haps();
        let nc = n_comp_haps as usize;
        for m in 0..self.n_markers {
            let mu = m as usize;
            let obs_allele = all_haps.allele(m, targ_hap);
            for j in 0..nc {
                if m == self.comp_hap_to_end[j] {
                    self.segment_index[j] += 1;
                    let si = self.segment_index[j];
                    self.comp_hap_to_hap[j] = self.comp_hap_hap[j].get(si);
                    self.comp_hap_to_end[j] = self.comp_hap_end[j].get(si);
                }
                let ref_hap = self.comp_hap_to_hap[j];
                haps[mu][j] = ref_hap;
                n_mismatches[mu][j] = u8::from(all_haps.allele(m, ref_hap) != obs_allele);
            }
        }
    }

    fn fill_q_with_random_haps(&mut self, hap: i32) {
        debug_assert!(self.q.is_empty());
        let ibs_haps = self.ibs_haps;
        let n_haps = ibs_haps.all_haps().n_haps();
        let n_states = (n_haps - 2).min(self.max_states);
        if n_states <= 0 {
            Utilities::exit("ERROR: there is only one sample");
        } else {
            let mut rand = Random::new(ibs_haps.phase_data().seed() + hap as i64);
            let sample = hap >> 1;
            let ibs_step = 0;
            let start_marker = 0;
            for j in 0..n_states {
                let mut h = rand.next_int_bound(n_haps);
                while (h >> 1) == sample {
                    h = rand.next_int_bound(n_haps);
                }
                let idx = self.q.size() as usize;
                self.comp_hap_hap[idx].add(h);
                self.q
                    .add(CompHapSegment::new(h, start_marker, ibs_step, j));
            }
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
    use std::rc::Rc;

    use super::super::{FixedPhaseData, PhaseData};

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_lfps_gt.vcf");
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
        let gts = [
            "0|0\t0|1\t1|1\t0|1\t0|0\t1|0",
            "0|1\t0|0\t1|0\t1|1\t0|1\t0|0",
            "1|1\t0|1\t0|0\t0|1\t1|0\t0|1",
            "0|0\t1|1\t0|1\t1|0\t0|0\t1|1",
            "0|1\t1|0\t1|1\t0|0\t0|1\t1|0",
            "1|0\t0|1\t0|0\t1|1\t1|0\t0|1",
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
    fn states_reference_valid_haps_and_consistent_mismatches() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let ibs = LowFreqPhaseIbs::new(&pd);
        let max_states = 4;
        let n_markers = pd.fpd().stage1_targ_gt().n_markers();
        let n_haps = pd.fpd().targ_gt().n_haps();
        let mut states = LowFreqPhaseStates::new(&ibs, max_states);

        let mut haps: Vec<Vec<i32>> = (0..n_markers)
            .map(|_| vec![0i32; max_states as usize])
            .collect();
        let mut mism: Vec<Vec<u8>> = (0..n_markers)
            .map(|_| vec![0u8; max_states as usize])
            .collect();

        let targ_hap = 2;
        let n_states = states.ibs_states(targ_hap, &mut haps, &mut mism);
        assert!(n_states >= 1 && n_states <= max_states);
        let all_haps = ibs.all_haps();
        for m in 0..n_markers as usize {
            for j in 0..n_states as usize {
                let rh = haps[m][j];
                assert!((0..n_haps).contains(&rh));
                // mismatch flag must agree with the phased reference alleles used to build it
                let obs = all_haps.allele(m as i32, targ_hap);
                let refa = all_haps.allele(m as i32, rh);
                assert_eq!(mism[m][j], u8::from(refa != obs));
            }
        }
    }
}
