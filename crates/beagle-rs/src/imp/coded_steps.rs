//! Port of `imp/CodedSteps.java` — partitions the imputation marker clusters into
//! non-overlapping steps (half-length first step, then `imp_step` cM apart), and for each step
//! indexes the distinct target allele sequences and maps every haplotype to its sequence index.
//! (Distinct from `phase::CodedSteps`.)

use crate::ints::IndexArray;

use super::ImpData;

/// Port of `imp/CodedSteps.java`.
pub struct CodedSteps<'a> {
    imp_data: &'a ImpData,
    step_starts: Vec<i32>,
    coded_steps: Vec<IndexArray>,
}

/// `Arrays.binarySearch(double[] a, int fromIndex, int toIndex, double key)`.
fn java_binary_search_f64(a: &[f64], from: usize, to: usize, key: f64) -> i32 {
    let mut low = from as i64;
    let mut high = to as i64 - 1;
    while low <= high {
        let mid = ((low + high) as u64 >> 1) as i64;
        let mid_val = a[mid as usize];
        if mid_val < key {
            low = mid + 1;
        } else if mid_val > key {
            high = mid - 1;
        } else {
            let mid_bits = mid_val.to_bits() as i64;
            let key_bits = key.to_bits() as i64;
            if mid_bits == key_bits {
                return mid as i32;
            } else if mid_bits < key_bits {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
    }
    -(low as i32 + 1)
}

fn next_index(pos: &[f64], start: usize, target_pos: f64) -> usize {
    let r = java_binary_search_f64(pos, start, pos.len(), target_pos);
    (if r < 0 { -r - 1 } else { r }) as usize
}

fn step_starts(imp_data: &ImpData) -> Vec<i32> {
    let pos = imp_data.pos_all();
    let step = imp_data.par().imp_step() as f64;
    let mut indices: Vec<i32> = Vec::with_capacity(pos.len() / 10);
    indices.push(0);
    let mut next_pos = pos[0] + step / 2.0; // make the first step half-length
    let mut index = next_index(&pos, 0, next_pos);
    while index < pos.len() {
        indices.push(index as i32);
        next_pos = pos[index] + step;
        index = next_index(&pos, index, next_pos);
    }
    indices
}

fn code_step(imp_data: &ImpData, starts: &[i32], start_index: usize) -> IndexArray {
    let n_ref_haps = imp_data.n_ref_haps();
    let n_haps = imp_data.n_haps();
    let mut hap_to_seq = vec![1i32; n_haps as usize];
    let start = starts[start_index];
    let end = if start_index + 1 < starts.len() {
        starts[start_index + 1]
    } else {
        imp_data.n_clusters()
    };
    let mut seq_cnt = 2; // seq 0 is reserved for sequences not found in the target
    for m in start..end {
        let h2s = imp_data.hap_to_seq(m);
        let coded_haps = h2s.int_array();
        let n_alleles = h2s.value_size();
        let mut seq_map = vec![0i32; (seq_cnt * n_alleles) as usize];
        seq_cnt = 1;
        for h in n_ref_haps..n_haps {
            let index = n_alleles * hap_to_seq[h as usize] + coded_haps.get(h);
            if seq_map[index as usize] == 0 {
                seq_map[index as usize] = seq_cnt;
                seq_cnt += 1;
            }
            hap_to_seq[h as usize] = seq_map[index as usize];
        }
        for h in 0..n_ref_haps {
            if hap_to_seq[h as usize] != 0 {
                let index = hap_to_seq[h as usize] * n_alleles + coded_haps.get(h);
                hap_to_seq[h as usize] = seq_map[index as usize];
            }
        }
    }
    IndexArray::from_ints(&hap_to_seq, seq_cnt)
}

impl<'a> CodedSteps<'a> {
    /// `new CodedSteps(ImpData impData)`.
    pub fn new(imp_data: &'a ImpData) -> Self {
        let step_starts = step_starts(imp_data);
        // Java parallelizes the per-step coding; the result order is preserved.
        let coded_steps = (0..step_starts.len())
            .map(|j| code_step(imp_data, &step_starts, j))
            .collect();
        CodedSteps {
            imp_data,
            step_starts,
            coded_steps,
        }
    }

    /// `impData()`.
    pub fn imp_data(&self) -> &ImpData {
        self.imp_data
    }
    /// `nSteps()`.
    pub fn n_steps(&self) -> i32 {
        self.step_starts.len() as i32
    }
    /// `stepStart(int step)`.
    pub fn step_start(&self, step: i32) -> i32 {
        self.step_starts[step as usize]
    }
    /// `get(int step)`.
    pub fn get(&self, step: i32) -> &IndexArray {
        &self.coded_steps[step as usize]
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
        let gt = std::env::temp_dir().join("beagle_rs_impcs_gt.vcf");
        let rf = std::env::temp_dir().join("beagle_rs_impcs_ref.vcf");
        std::fs::write(&gt, b"x").unwrap();
        std::fs::write(&rf, b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            format!("ref={}", rf.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
        ]);
        let hr = header(3);
        let ht = header(2);
        let mp = MarkerParser::new(true, true, true, true);
        let ref_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|1\t0|1\t1|0",
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
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|1",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|0\t1|0",
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
    fn codes_steps_over_clusters() {
        let imp = imp_data();
        let cs = CodedSteps::new(&imp);
        assert!(cs.n_steps() >= 1);
        assert_eq!(cs.step_start(0), 0);
        let n_haps = imp.n_haps();
        for step in 0..cs.n_steps() {
            let ia = cs.get(step);
            assert_eq!(ia.int_array().size(), n_haps);
            let vs = ia.value_size();
            for h in 0..n_haps {
                let seq = ia.int_array().get(h);
                assert!(seq >= 0 && seq < vs);
            }
        }
    }
}
