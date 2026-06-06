//! Port of `vcf/RestrictedGT.java` — a `GT` wrapper that restricts an underlying `GT` to a
//! subset of its markers via an inclusion map.
//!
//! Note: Java's `toString()` is `"class vcf.RestrictedGT : " + gt.toString()`. `GT`
//! toString is debug-only (and non-deterministic for `RefGT`, whose toString is the default
//! `Object` form) and never part of byte-exact output, so the Rust `Display` omits the
//! inner record dump.

use std::fmt;
use std::rc::Rc;

use super::{Marker, Markers, Samples, GT};

/// Port of `vcf/RestrictedGT.java`.
pub struct RestrictedGT {
    gt: Rc<dyn GT>,
    restricted_markers: Markers,
    inclusion_map: Vec<i32>,
}

impl RestrictedGT {
    /// `new RestrictedGT(GT gt, Markers markers, int[] indices)`.
    pub fn new(gt: Rc<dyn GT>, markers: Markers, indices: &[i32]) -> Self {
        for (j, &idx) in indices.iter().enumerate() {
            assert!(!(j > 0 && idx <= indices[j - 1]), "{idx}");
            assert!(
                gt.marker(idx) == markers.marker(j as i32),
                "{}",
                markers.marker(j as i32)
            );
        }
        RestrictedGT {
            gt,
            restricted_markers: markers,
            inclusion_map: indices.to_vec(),
        }
    }
}

impl GT for RestrictedGT {
    fn is_reversed(&self) -> bool {
        false
    }

    fn n_markers(&self) -> i32 {
        self.restricted_markers.size()
    }

    fn marker(&self, marker: i32) -> &Marker {
        self.restricted_markers.marker(marker)
    }

    fn markers(&self) -> &Markers {
        &self.restricted_markers
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
        self.gt.allele(self.inclusion_map[marker as usize], hap)
    }

    fn restrict(self: Rc<Self>, markers: &Markers, indices: &[i32]) -> Rc<dyn GT> {
        Rc::new(RestrictedGT::new(self, markers.clone(), indices))
    }

    fn restrict_range(self: Rc<Self>, start: i32, end: i32) -> Rc<dyn GT> {
        let restrict_markers = self.restricted_markers.restrict(start, end);
        let indices: Vec<i32> = (start..end).collect();
        Rc::new(RestrictedGT::new(self, restrict_markers, &indices))
    }
}

impl fmt::Display for RestrictedGT {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Java: "class vcf.RestrictedGT : " + gt.toString(); inner omitted (debug-only).
        f.write_str("class vcf.RestrictedGT")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::GTRec;
    use crate::vcf::{BasicGT, BasicGTRec, MarkerParser, VcfHeader, VcfRecGTParser, HEADER_PREFIX};

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

    fn base_gt() -> Rc<BasicGT> {
        let h = header(1);
        Rc::new(BasicGT::new(vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1"),
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|0"),
            rec(&h, "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|0"),
        ]))
    }

    #[test]
    fn restricts_to_subset_via_inclusion_map() {
        let gt = base_gt();
        let markers = gt.markers().restrict_indices(&[0, 2]);
        let r = RestrictedGT::new(gt.clone(), markers, &[0, 2]);
        assert_eq!(r.n_markers(), 2);
        assert_eq!(r.marker(0).pos(), 100);
        assert_eq!(r.marker(1).pos(), 300);
        // allele(restrictedMarker, hap) maps through the inclusion map to the original
        assert_eq!(r.allele(0, 0), gt.allele(0, 0)); // marker 0 -> orig 0
        assert_eq!(r.allele(1, 1), gt.allele(2, 1)); // marker 1 -> orig 2
        assert_eq!(r.n_haps(), gt.n_haps());
        assert!(r.is_phased());
    }

    #[test]
    fn nested_restrict_range() {
        let gt = base_gt();
        let markers = gt.markers().clone();
        let r: Rc<dyn GT> = Rc::new(RestrictedGT::new(gt, markers, &[0, 1, 2]));
        let sub = r.restrict_range(1, 3);
        assert_eq!(sub.n_markers(), 2);
        assert_eq!(sub.marker(0).pos(), 200);
        assert_eq!(sub.marker(1).pos(), 300);
    }

    #[test]
    #[should_panic(expected = "0")]
    fn rejects_non_increasing_indices() {
        let gt = base_gt();
        let markers = gt.markers().clone();
        let _ = RestrictedGT::new(gt, markers, &[0, 0, 2]);
    }
}
