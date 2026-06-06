//! Port of `imp/ImpLS.java` — computes the HMM state probabilities at the genotyped markers for
//! every target haplotype (one `ImpLSBaum::impute` per haplotype).

use super::{ImpData, ImpIbs, ImpLSBaum, StateProbs};

/// Port of `imp/ImpLS.java`.
pub struct ImpLS;

impl ImpLS {
    /// `stateProbs(ImpData impData)` — HMM state probabilities per target haplotype.
    ///
    /// Java distributes haplotypes across `nthreads` workers sharing one `ImpIbs`; each
    /// haplotype's result is deterministic (RNG seeded by the hap index), so the
    /// single-threaded port computes them in order with one `ImpLSBaum`.
    pub fn state_probs(imp_data: &ImpData) -> Vec<Box<dyn StateProbs>> {
        let n_targ_haps = imp_data.targ_gt().n_haps();
        let ibs_haps = ImpIbs::new(imp_data);
        let mut baum = ImpLSBaum::new(imp_data, &ibs_haps);
        let mut st_probs: Vec<Box<dyn StateProbs>> = Vec::with_capacity(n_targ_haps as usize);
        for hap in 0..n_targ_haps {
            st_probs.push(Box::new(baum.impute(hap)));
        }
        st_probs
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
        let gt = std::env::temp_dir().join("beagle_rs_impls_gt.vcf");
        let rf = std::env::temp_dir().join("beagle_rs_impls_ref.vcf");
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
    fn state_probs_per_target_hap() {
        let imp = imp_data();
        let n_targ_haps = imp.targ_gt().n_haps();
        let n_clusters = imp.n_clusters();
        let probs = ImpLS::state_probs(&imp);
        assert_eq!(probs.len() as i32, n_targ_haps);
        for (h, sp) in probs.iter().enumerate() {
            assert_eq!(sp.targ_hap(), h as i32);
            assert_eq!(sp.n_targ_markers(), n_clusters);
        }
    }
}
