//! Port of `imp/ImpIbs.java` — for each step and target haplotype, finds reference haplotypes
//! sharing a long IBS segment, by recursively partitioning haplotypes by their coded allele
//! sequence over successive steps until a partition has few enough reference haplotypes, then
//! taking (and randomly padding/truncating to) `nHapsPerStep` IBS reference haplotypes.

use std::rc::Rc;

use crate::blbutil::Utilities;
use crate::ints::{IndexArray, IntList};
use crate::jdk::Random;

use super::{CodedSteps, ImpData};

/// Port of `imp/ImpIbs.java`.
pub struct ImpIbs<'a> {
    imp_data: &'a ImpData,
    n_ref_haps: i32,
    seed: i64,
    n_steps: i32,
    n_haps_per_step: i32,
    coded_steps: CodedSteps<'a>,
    ibs_haps: Vec<Vec<Rc<Vec<i32>>>>, // [step][targ hap] -> shared IBS set
}

fn ins_pt(il: &IntList, n_ref_haps: i32) -> i32 {
    let index = il.binary_search(n_ref_haps);
    if index >= 0 {
        index
    } else {
        -index - 1
    }
}

fn random_subset(list: &IntList, size: i32, rand: &mut Random) -> Vec<i32> {
    let size = size.min(list.size());
    let mut ia = list.to_array();
    for j in 0..size as usize {
        let x = rand.next_int_bound(ia.len() as i32 - j as i32) as usize;
        ia.swap(j, j + x);
    }
    ia.truncate(size as usize);
    ia
}

impl<'a> ImpIbs<'a> {
    /// `new ImpIbs(ImpData impData)`.
    pub fn new(imp_data: &'a ImpData) -> Self {
        let par = imp_data.par();
        let seed = par.seed();
        let n_ref_haps = imp_data.n_ref_haps();
        let n_steps = par.imp_nsteps();
        // Math.round(imp_segment/imp_step) = floor(x + 0.5)
        let n_steps_per_segment = (par.imp_segment() / par.imp_step() + 0.5).floor() as i32;
        let n_haps_per_step = par.imp_states() / n_steps_per_segment;
        let coded_steps = CodedSteps::new(imp_data);
        let this = ImpIbs {
            imp_data,
            n_ref_haps,
            seed,
            n_steps,
            n_haps_per_step,
            coded_steps,
            ibs_haps: Vec::new(),
        };
        let n = this.coded_steps.n_steps();
        // Java parallelizes per step; the result order is preserved.
        let ibs_haps = (0..n).map(|j| this.get_ibs_haps(j)).collect();
        ImpIbs { ibs_haps, ..this }
    }

    fn get_ibs_haps(&self, index: i32) -> Vec<Rc<Vec<i32>>> {
        let n_targ_haps = self.imp_data.n_targ_haps();
        let mut results: Vec<Option<Rc<Vec<i32>>>> = vec![None; n_targ_haps as usize];
        let n_steps_to_merge = self.n_steps.min(self.coded_steps.n_steps() - index);
        let children = self.init_partition(self.coded_steps.get(index));
        let mut next_parents: Vec<IntList> = Vec::with_capacity(children.len());
        self.init_update_results(children, &mut next_parents, &mut results);
        for i in 1..n_steps_to_merge {
            let parents = std::mem::take(&mut next_parents);
            let init_capacity = n_targ_haps.min(2 * parents.len() as i32).max(0);
            next_parents = Vec::with_capacity(init_capacity as usize);
            let coded_step = self.coded_steps.get(index + i);
            for parent in &parents {
                let children = self.partition(parent, coded_step);
                self.update_results(parent, children, &mut next_parents, &mut results);
            }
        }
        self.final_update_results(next_parents, &mut results);
        results.into_iter().map(|r| r.expect("ibs set")).collect()
    }

    fn init_partition(&self, coded_step: &IndexArray) -> Vec<IntList> {
        let mut list: Vec<Option<usize>> = vec![None; coded_step.value_size() as usize];
        let hap2_seq = coded_step.int_array();
        let n_haps = hap2_seq.size();
        let mut children: Vec<IntList> = Vec::new();
        for h in self.n_ref_haps..n_haps {
            let seq = hap2_seq.get(h) as usize;
            if list[seq].is_none() {
                list[seq] = Some(children.len());
                children.push(IntList::new());
            }
        }
        for h in 0..n_haps {
            let seq = hap2_seq.get(h) as usize;
            if let Some(idx) = list[seq] {
                children[idx].add(h);
            }
        }
        children
    }

    fn partition(&self, parent: &IntList, coded_step: &IndexArray) -> Vec<IntList> {
        let mut list: Vec<Option<usize>> = vec![None; coded_step.value_size() as usize];
        let hap2_seq = coded_step.int_array();
        let n_parent_haps = parent.size();
        let mut children: Vec<IntList> = Vec::new();
        let targ_start = ins_pt(parent, self.n_ref_haps);
        for k in targ_start..n_parent_haps {
            let hap = parent.get(k);
            let seq = hap2_seq.get(hap) as usize;
            if list[seq].is_none() {
                list[seq] = Some(children.len());
                children.push(IntList::new());
            }
        }
        for k in 0..n_parent_haps {
            let hap = parent.get(k);
            let seq = hap2_seq.get(hap) as usize;
            if let Some(idx) = list[seq] {
                children[idx].add(hap);
            }
        }
        children
    }

    fn init_update_results(
        &self,
        children: Vec<IntList>,
        next_parents: &mut Vec<IntList>,
        results: &mut [Option<Rc<Vec<i32>>>],
    ) {
        for hap_list in children {
            let n_ref = ins_pt(&hap_list, self.n_ref_haps);
            if n_ref <= self.n_haps_per_step {
                let ibs = Rc::new(hap_list.copy_of(n_ref));
                self.set_result(&hap_list, n_ref, ibs, results);
            } else {
                next_parents.push(hap_list);
            }
        }
    }

    fn update_results(
        &self,
        parent: &IntList,
        children: Vec<IntList>,
        next_ibs: &mut Vec<IntList>,
        results: &mut [Option<Rc<Vec<i32>>>],
    ) {
        for child in children {
            let n_child_ref = ins_pt(&child, self.n_ref_haps);
            if n_child_ref <= self.n_haps_per_step {
                let ibs_list = self.ibs_haps_set(parent, &child, n_child_ref);
                self.set_result(&child, n_child_ref, Rc::new(ibs_list), results);
            } else {
                next_ibs.push(child);
            }
        }
    }

    fn ibs_haps_set(&self, parent: &IntList, child: &IntList, n_child_ref: i32) -> Vec<i32> {
        let mut combined = IntList::with_capacity(self.n_haps_per_step);
        for j in 0..n_child_ref {
            combined.add(child.get(j));
        }
        let size = self.n_haps_per_step - n_child_ref;
        let mut rand = Random::new(self.seed + parent.get(0) as i64);
        let uniq_to_parent = self.uniq_to_parent(parent, child, n_child_ref);
        let rand_subset = random_subset(&uniq_to_parent, size, &mut rand);
        for i in rand_subset {
            combined.add(i);
        }
        let mut ia = combined.to_array();
        ia.sort_unstable();
        ia
    }

    fn uniq_to_parent(&self, parent: &IntList, child: &IntList, n_child_ref: i32) -> IntList {
        let n_child_ref_m1 = n_child_ref - 1;
        let n_parent_ref = ins_pt(parent, self.n_ref_haps);
        let mut uniq = IntList::with_capacity(parent.size());
        let mut c = 0;
        let mut c_val = child.get(c);
        for p in 0..n_parent_ref {
            let p_val = parent.get(p);
            while c_val < p_val && c < n_child_ref_m1 {
                c += 1;
                c_val = child.get(c);
            }
            if p_val != c_val {
                uniq.add(p_val);
            }
        }
        uniq
    }

    fn final_update_results(&self, children: Vec<IntList>, results: &mut [Option<Rc<Vec<i32>>>]) {
        for child in children {
            let n_ref = ins_pt(&child, self.n_ref_haps);
            let mut ibs_list = child.copy_of(n_ref);
            if self.n_haps_per_step < ibs_list.len() as i32 {
                let mut rand = Random::new(self.seed + child.get(0) as i64);
                Utilities::shuffle(&mut ibs_list, &mut rand);
                ibs_list.truncate(self.n_haps_per_step as usize);
                ibs_list.sort_unstable();
            }
            self.set_result(&child, n_ref, Rc::new(ibs_list), results);
        }
    }

    fn set_result(
        &self,
        child: &IntList,
        first_targ_index: i32,
        ibs_haps: Rc<Vec<i32>>,
        results: &mut [Option<Rc<Vec<i32>>>],
    ) {
        for j in first_targ_index..child.size() {
            results[(child.get(j) - self.n_ref_haps) as usize] = Some(ibs_haps.clone());
        }
    }

    /// `ibsHaps(int hap, int step)` — reference haplotype indices IBS with target `hap` over the
    /// interval starting at `step`.
    pub fn ibs_haps(&self, hap: i32, step: i32) -> Vec<i32> {
        (*self.ibs_haps[step as usize][hap as usize]).clone()
    }

    /// `impData()`.
    pub fn imp_data(&self) -> &ImpData {
        self.imp_data
    }

    /// `codedSteps()`.
    pub fn coded_steps(&self) -> &CodedSteps<'a> {
        &self.coded_steps
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let gt = std::env::temp_dir().join("beagle_rs_impibs_gt.vcf");
        let rf = std::env::temp_dir().join("beagle_rs_impibs_ref.vcf");
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
    fn finds_ibs_ref_haps_per_target_hap() {
        let imp = imp_data();
        let n_ref_haps = imp.n_ref_haps();
        let n_targ_haps = imp.n_targ_haps();
        let ibs = ImpIbs::new(&imp);
        let n_steps = ibs.coded_steps().n_steps();
        for step in 0..n_steps {
            for hap in 0..n_targ_haps {
                let haps = ibs.ibs_haps(hap, step);
                // every returned index is a reference haplotype, sorted, within capacity
                for &h in &haps {
                    assert!(h >= 0 && h < n_ref_haps, "step {step} hap {hap} -> {h}");
                }
                assert!(haps.windows(2).all(|w| w[0] <= w[1]), "not sorted");
            }
        }
    }
}
