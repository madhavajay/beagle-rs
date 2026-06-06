//! Port of `imp/ImpStates.java` — builds up to `maxStates` pseudo-reference haplotypes (mosaics
//! of reference-haplotype segments) for a target haplotype from the `ImpIbs` IBS segments, using
//! a `lastIbsStep` priority queue, then emits per-cluster state reference haps and allele-match
//! flags.

use crate::beagleutil::CompHapSegment;
use crate::ints::{IntIntMap, IntList};
use crate::jdk::{PriorityQueue, Random};

use super::ImpIbs;

const NIL: i32 = -103;

/// Port of `imp/ImpStates.java`.
pub struct ImpStates<'a> {
    ibs_haps: &'a ImpIbs<'a>,
    n_clusters: i32,
    max_states: i32,
    hap_to_last_ibs_step: IntIntMap,
    q: PriorityQueue<CompHapSegment>,
    comp_hap_hap: Vec<IntList>,
    comp_hap_end: Vec<IntList>,
    comp_hap_to_list_index: Vec<i32>,
    comp_hap_to_hap: Vec<i32>,
    comp_hap_to_end: Vec<i32>,
}

impl<'a> ImpStates<'a> {
    /// `new ImpStates(ImpIbs ibsHaps)`.
    pub fn new(ibs_haps: &'a ImpIbs<'a>) -> Self {
        let imp_data = ibs_haps.imp_data();
        let n_clusters = imp_data.n_clusters();
        let max_states = imp_data.par().imp_states();
        let ms = max_states as usize;
        ImpStates {
            ibs_haps,
            n_clusters,
            max_states,
            hap_to_last_ibs_step: IntIntMap::new(max_states),
            q: PriorityQueue::new(ms),
            comp_hap_hap: (0..max_states).map(|_| IntList::new()).collect(),
            comp_hap_end: (0..max_states).map(|_| IntList::new()).collect(),
            comp_hap_to_list_index: vec![0; ms],
            comp_hap_to_hap: vec![0; ms],
            comp_hap_to_end: vec![0; ms],
        }
    }

    /// `maxStates()`.
    pub fn max_states(&self) -> i32 {
        self.max_states
    }

    /// `ibsStates(int targHap, int[][] haps, boolean[][] alMatch)` — fills per-cluster state
    /// reference haps and allele-match flags, returning the number of states.
    pub fn ibs_states(
        &mut self,
        targ_hap: i32,
        haps: &mut [Vec<i32>],
        al_match: &mut [Vec<bool>],
    ) -> i32 {
        self.initialize_fields();
        let ibs_haps = self.ibs_haps;
        let n = ibs_haps.coded_steps().n_steps();
        for j in 0..n {
            let ibs = ibs_haps.ibs_haps(targ_hap, j);
            for hap in ibs {
                self.update_fields(hap, j);
            }
        }
        if self.q.is_empty() {
            self.fill_q_with_random_haps(targ_hap);
        }
        self.copy_data(targ_hap, haps, al_match)
    }

    fn initialize_fields(&mut self) {
        self.hap_to_last_ibs_step.clear();
        for j in 0..self.q.size() as usize {
            self.comp_hap_hap[j].clear();
            self.comp_hap_end[j].clear();
        }
        self.q.clear();
    }

    fn update_fields(&mut self, hap: i32, step: i32) {
        if self.hap_to_last_ibs_step.get(hap, NIL) == NIL {
            // hap not currently in q
            self.update_head_of_q();
            if self.q.size() == self.max_states {
                let mut head = self.q.poll().expect("nonempty");
                let start_marker = self
                    .ibs_haps
                    .coded_steps()
                    .step_start((head.last_ibs_step() + step) >> 1);
                self.hap_to_last_ibs_step.remove(head.hap());
                let idx = head.comp_hap_index() as usize;
                self.comp_hap_hap[idx].add(hap); // hap of new segment
                self.comp_hap_end[idx].add(start_marker); // end of previous segment
                head.update_segment(hap, start_marker, step);
                self.q.offer(head);
            } else {
                let comp_hap_index = self.q.size();
                let start_marker = 0;
                self.comp_hap_hap[comp_hap_index as usize].add(hap);
                self.q
                    .offer(CompHapSegment::new(hap, start_marker, step, comp_hap_index));
            }
        }
        self.hap_to_last_ibs_step.put(hap, step);
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

    fn copy_data(
        &mut self,
        targ_hap: i32,
        hap_indices: &mut [Vec<i32>],
        al_match: &mut [Vec<bool>],
    ) -> i32 {
        let n_comp_haps = self.q.size();
        let imp_data = self.ibs_haps.imp_data();
        let shifted_targ_hap = imp_data.n_ref_haps() + targ_hap;
        self.initialize_copy(n_comp_haps);
        for m in 0..self.n_clusters {
            let targ_allele = imp_data.allele(m, shifted_targ_hap);
            for j in 0..n_comp_haps as usize {
                if m == self.comp_hap_to_end[j] {
                    self.comp_hap_to_list_index[j] += 1;
                    let li = self.comp_hap_to_list_index[j];
                    self.comp_hap_to_hap[j] = self.comp_hap_hap[j].get(li);
                    self.comp_hap_to_end[j] = self.comp_hap_end[j].get(li);
                }
                hap_indices[m as usize][j] = self.comp_hap_to_hap[j];
                al_match[m as usize][j] =
                    imp_data.allele(m, self.comp_hap_to_hap[j]) == targ_allele;
            }
        }
        n_comp_haps
    }

    fn initialize_copy(&mut self, n_slots: i32) {
        for j in 0..n_slots as usize {
            self.comp_hap_end[j].add(self.n_clusters); // add missing end of last segment
            self.comp_hap_to_list_index[j] = 0;
            self.comp_hap_to_hap[j] = self.comp_hap_hap[j].get(0);
            self.comp_hap_to_end[j] = self.comp_hap_end[j].get(0);
        }
    }

    fn fill_q_with_random_haps(&mut self, hap: i32) {
        debug_assert!(self.q.is_empty());
        let n_ref_haps = self.ibs_haps.imp_data().n_ref_haps();
        let n_states = n_ref_haps.min(self.max_states);
        let mut rand = Random::new(hap as i64);
        let ibs_step = 0;
        let start_marker = 0;
        for j in 0..n_states {
            let mut h = rand.next_int_bound(n_ref_haps);
            while h == hap {
                h = rand.next_int_bound(n_ref_haps);
            }
            self.comp_hap_hap[j as usize].add(h);
            self.q
                .add(CompHapSegment::new(h, start_marker, ibs_step, j));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::imp::ImpData;
    use crate::main_pkg::Par;
    use crate::vcf::{
        allele_ref_gt_rec_from_parser, BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices,
        MarkerParser, PositionMap, RefGT, RefGTRec, VcfHeader, VcfRecGTParser, Window, GT,
        HEADER_PREFIX,
    };
    use std::rc::Rc;

    fn header(n: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n {
            hdr.push_str(&format!("\tS{s}"));
        }
        VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &vec![true; n],
        )
    }

    fn imp_data() -> ImpData {
        let gt = std::env::temp_dir().join("beagle_rs_impstates_gt.vcf");
        let rf = std::env::temp_dir().join("beagle_rs_impstates_ref.vcf");
        std::fs::write(&gt, b"x").unwrap();
        std::fs::write(&rf, b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            format!("ref={}", rf.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
        ]);
        let hr = header(6);
        let ht = header(3);
        let mp = MarkerParser::new(true, true, true, true);
        let ref_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|1\t1|0\t0|0",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\t1|1\t0|1\t1|0",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1\t0|0\t1|0\t0|1",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|1\t0|1\t1|0\t1|1\t0|0\t1|0",
        ];
        let ref_recs: Vec<Rc<dyn RefGTRec>> = ref_lines
            .iter()
            .map(|l| {
                Rc::from(allele_ref_gt_rec_from_parser(&VcfRecGTParser::new(
                    &hr, l, &mp,
                )))
            })
            .collect();
        let ref_gt = RefGT::from_recs(ref_recs);
        let targ_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|0\t1|0\t0|1",
        ];
        let targ_recs: Vec<Rc<dyn GTRec>> = targ_lines
            .iter()
            .map(|l| {
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(&ht, l, &mp))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(targ_recs);
        let indices = MarkerIndices::from_in_targ(&[true, false, false, true], 0, 4);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        let window = Window::new(gen_map.clone(), 1, true, indices, Some(ref_gt), targ);
        let phased: Rc<dyn GT> = Rc::new(window.targ_gt().clone());
        ImpData::new(&par, &window, phased, gen_map.as_ref())
    }

    #[test]
    fn builds_states_and_allele_match() {
        let imp = imp_data();
        let n_ref_haps = imp.n_ref_haps();
        let n_clusters = imp.n_clusters();
        let ibs = ImpIbs::new(&imp);
        let mut states = ImpStates::new(&ibs);
        let max_states = states.max_states();
        let mut haps: Vec<Vec<i32>> = (0..n_clusters)
            .map(|_| vec![0i32; max_states as usize])
            .collect();
        let mut al_match: Vec<Vec<bool>> = (0..n_clusters)
            .map(|_| vec![false; max_states as usize])
            .collect();
        let n_states = states.ibs_states(0, &mut haps, &mut al_match);
        assert!(n_states >= 1 && n_states <= max_states);
        // every state hap is a valid reference hap; al_match is consistent with the alleles
        for m in 0..n_clusters as usize {
            for j in 0..n_states as usize {
                let rh = haps[m][j];
                assert!(rh >= 0 && rh < n_ref_haps);
                let expected = imp.allele(m as i32, rh) == imp.allele(m as i32, n_ref_haps);
                assert_eq!(al_match[m][j], expected);
            }
        }
    }
}
