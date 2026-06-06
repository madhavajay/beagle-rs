//! Port of `phase/CodedSteps.java` — partitions the phased haplotypes into the
//! non-overlapping `Steps` intervals, indexes the distinct allele sequence in each interval,
//! and stores, per step, a haplotype→sequence-index map (an `IndexArray`).

use crate::ints::{IndexArray, IntIntMap};
use crate::vcf::{Samples, Steps, XRefGT, GT};

use super::EstPhase;

/// Port of `phase/CodedSteps.java`.
pub struct CodedSteps {
    targ_samples: Samples,
    ref_haps: Option<XRefGT>,
    all_haps: XRefGT,
    steps: Steps,
    coded_steps: Vec<IndexArray>,
}

impl CodedSteps {
    /// `new CodedSteps(EstPhase estPhase)`.
    pub fn new(est_phase: &EstPhase) -> Self {
        let fpd = est_phase.fpd();
        let targ_haps = est_phase.phased_haps();
        let targ_samples = targ_haps.samples().clone();
        let ref_haps = fpd.stage1_xref_gt().cloned();
        let all_haps = match &ref_haps {
            Some(rh) => XRefGT::combine(&targ_haps, rh),
            None => targ_haps,
        };
        let steps = fpd.stage1_steps().clone();
        let coded_steps = coded_steps(&all_haps, &steps);
        CodedSteps {
            targ_samples,
            ref_haps,
            all_haps,
            steps,
            coded_steps,
        }
    }

    /// `get(int step)` — the haplotype→sequence-index map for `step`.
    pub fn get(&self, step: i32) -> &IndexArray {
        &self.coded_steps[step as usize]
    }

    /// `targSamples()`.
    pub fn targ_samples(&self) -> &Samples {
        &self.targ_samples
    }

    /// `refHaps()`.
    pub fn ref_haps(&self) -> Option<&XRefGT> {
        self.ref_haps.as_ref()
    }

    /// `allHaps()` — phased target+reference haplotypes (target haps first).
    pub fn all_haps(&self) -> &XRefGT {
        &self.all_haps
    }

    /// `steps()`.
    pub fn steps(&self) -> &Steps {
        &self.steps
    }
}

/// `codedSteps(XRefGT gt, Steps steps, int nThreads)` — Java batches the steps across threads
/// and concatenates the per-batch results in step order; each step is indexed independently,
/// so the sequential build over `[0, steps.size())` produces the identical array.
fn coded_steps(gt: &XRefGT, steps: &Steps) -> Vec<IndexArray> {
    let sentinel = -1;
    let n_haps = gt.n_haps();
    let n_steps = steps.size();
    (0..n_steps)
        .map(|j| {
            let m_start = steps.start(j);
            let m_end = steps.end(j);
            let mut hap_to_seq = vec![0i32; n_haps as usize];
            let mut seq_cnt = 0i32;
            let mut seq_map = IntIntMap::new(8);
            for h in 0..n_haps {
                let key = gt.hash(h, m_start, m_end);
                let mut seq_index = seq_map.get(key, sentinel);
                if seq_index == sentinel {
                    seq_index = seq_cnt;
                    seq_cnt += 1;
                    seq_map.put(key, seq_index);
                }
                hap_to_seq[h as usize] = seq_index;
            }
            IndexArray::from_ints(&hap_to_seq, seq_cnt)
        })
        .collect()
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
        let gt = std::env::temp_dir().join("beagle_rs_cs_gt.vcf");
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
    fn indexes_haplotype_sequences_per_step() {
        let est = EstPhase::new(Rc::new(fpd()), 99999);
        let cs = CodedSteps::new(&est);
        assert_eq!(cs.targ_samples().size(), 3);
        assert!(cs.ref_haps().is_none());
        assert_eq!(cs.all_haps().n_haps(), 6); // 3 diploid samples, no ref

        // Each step maps all 6 haplotypes to a sequence index in [0, valueSize).
        for step in 0..cs.steps().size() {
            let ia = cs.get(step);
            assert_eq!(ia.int_array().size(), 6);
            let vs = ia.value_size();
            assert!(vs >= 1);
            for h in 0..6 {
                let seq = ia.int_array().get(h);
                assert!(seq >= 0 && seq < vs);
            }
            // S0 (haps 0,1) is all-ref "0|0" and S2 (haps 4,5) is all-alt "1|1": within a
            // step these are distinct sequences, so at least 2 indices are used.
            assert!(vs >= 2);
        }
    }
}
