//! Port of `vcf/SplicedGT.java` — genotypes obtained by replacing the initial `overlap`
//! markers of one `GT` with phased genotypes from another `GT`.
//!
//! Note: the Java constructor's javadoc mentions an `isPhased` precondition but the code
//! does not check it; the Rust port matches (no `isPhased` check).

use std::fmt;
use std::rc::Rc;

use super::{Marker, Markers, RestrictedGT, Samples, GT};

/// Port of `vcf/SplicedGT.java`.
pub struct SplicedGT {
    overlap: i32,
    phased_gt: Rc<dyn GT>,
    gt: Rc<dyn GT>,
}

impl SplicedGT {
    /// `new SplicedGT(GT phasedOverlap, GT gt)`.
    pub fn new(phased_overlap: Rc<dyn GT>, gt: Rc<dyn GT>) -> Self {
        assert!(
            phased_overlap.n_markers() < gt.n_markers(),
            "inconsistent markers"
        );
        for j in 0..phased_overlap.n_markers() {
            assert!(
                phased_overlap.marker(j) == gt.marker(j),
                "inconsistent markers"
            );
        }
        assert!(
            phased_overlap.samples() == gt.samples(),
            "inconsistent samples"
        );
        SplicedGT {
            overlap: phased_overlap.n_markers(),
            phased_gt: phased_overlap,
            gt,
        }
    }
}

impl GT for SplicedGT {
    fn is_reversed(&self) -> bool {
        false
    }

    fn marker(&self, marker: i32) -> &Marker {
        self.gt.marker(marker)
    }

    fn markers(&self) -> &Markers {
        self.gt.markers()
    }

    fn n_markers(&self) -> i32 {
        self.gt.n_markers()
    }

    fn n_haps(&self) -> i32 {
        self.gt.n_haps()
    }

    fn n_samples(&self) -> i32 {
        self.gt.n_samples()
    }

    fn samples(&self) -> &Samples {
        self.gt.samples()
    }

    fn is_phased(&self) -> bool {
        self.gt.is_phased()
    }

    fn allele(&self, marker: i32, hap: i32) -> i32 {
        if marker < self.overlap {
            self.phased_gt.allele(marker, hap)
        } else {
            self.gt.allele(marker, hap)
        }
    }

    fn restrict(self: Rc<Self>, markers: &Markers, indices: &[i32]) -> Rc<dyn GT> {
        Rc::new(RestrictedGT::new(self, markers.clone(), indices))
    }

    fn restrict_range(self: Rc<Self>, start: i32, end: i32) -> Rc<dyn GT> {
        if start >= self.overlap {
            self.gt.clone().restrict_range(start, end)
        } else {
            let restrict_phased = self.phased_gt.clone().restrict_range(start, self.overlap);
            let restrict_gt = self.gt.clone().restrict_range(start, end);
            Rc::new(SplicedGT::new(restrict_phased, restrict_gt))
        }
    }
}

impl fmt::Display for SplicedGT {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SplicedGL: nSamples={}", self.n_samples())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, MarkerParser, VcfHeader, VcfRecGTParser, HEADER_PREFIX,
    };

    fn header(n_dip: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n_dip {
            hdr.push_str(&format!("\tS{s}"));
        }
        let lines = vec!["##fileformat=VCFv4.2".to_string(), hdr];
        VcfHeader::new_accept_all("src", &lines, &vec![true; n_dip])
    }

    fn rec(h: &VcfHeader, line: &str) -> Rc<dyn GTRec> {
        let p = VcfRecGTParser::new(h, line, &MarkerParser::new(true, true, true, true));
        Rc::new(BasicGTRec::from_parser(&p))
    }

    #[test]
    fn splices_initial_markers_from_phased_gt() {
        let h = header(1);
        // Same markers (chrom/pos/alleles) but different genotype values, so allele()
        // distinguishes the two sources.
        let phased: Rc<dyn GT> = Rc::new(BasicGT::new(vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t1|1"),
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|1"),
        ]));
        let full: Rc<dyn GT> = Rc::new(BasicGT::new(vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0"),
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0"),
            rec(&h, "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|0"),
        ]));
        let s = SplicedGT::new(phased, full);
        assert_eq!(s.n_markers(), 3);
        // markers 0,1 come from phased (allele 1); marker 2 from gt (allele 0)
        assert_eq!(s.allele(0, 0), 1);
        assert_eq!(s.allele(1, 1), 1);
        assert_eq!(s.allele(2, 0), 0);
        assert_eq!(s.to_string(), "SplicedGL: nSamples=1");
    }

    #[test]
    fn restrict_range_after_overlap_returns_inner_gt() {
        let h = header(1);
        let phased: Rc<dyn GT> = Rc::new(BasicGT::new(vec![rec(
            &h,
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t1|1",
        )]));
        let full: Rc<dyn GT> = Rc::new(BasicGT::new(vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0"),
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1"),
            rec(&h, "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|0"),
        ]));
        let s: Rc<dyn GT> = Rc::new(SplicedGT::new(phased, full));
        // overlap == 1; restrict [1,3) starts at overlap -> delegates to inner gt
        let sub = s.restrict_range(1, 3);
        assert_eq!(sub.n_markers(), 2);
        assert_eq!(sub.allele(0, 1), 1); // gt marker 1 hap 1
    }

    #[test]
    fn restrict_range_spanning_overlap_returns_spliced() {
        let h = header(1);
        let phased: Rc<dyn GT> = Rc::new(BasicGT::new(vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t1|1"),
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|1"),
        ]));
        let full: Rc<dyn GT> = Rc::new(BasicGT::new(vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0"),
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0"),
            rec(&h, "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|0"),
        ]));
        let s: Rc<dyn GT> = Rc::new(SplicedGT::new(phased, full));
        // overlap == 2; restrict [0,3) spans the overlap -> still spliced
        let sub = s.restrict_range(0, 3);
        assert_eq!(sub.n_markers(), 3);
        assert_eq!(sub.allele(0, 0), 1); // from phased
        assert_eq!(sub.allele(2, 0), 0); // from gt
    }
}
