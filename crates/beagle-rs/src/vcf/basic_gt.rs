//! Port of `vcf/BasicGT.java` — genotypes for a list of markers and samples, backed by
//! one `GTRec` per marker. Records are shared (via `Rc`) when a `BasicGT` is restricted.

use crate::blbutil::consts;
use std::fmt;
use std::rc::Rc;

use super::{to_vcf_rec, GTRec, Marker, Markers, Samples, GT};

/// Port of `vcf/BasicGT.java`.
pub struct BasicGT {
    samples: Samples,
    markers: Markers,
    recs: Vec<Rc<dyn GTRec>>,
    is_ref_data: bool,
}

/// `BasicGT.genotype(int a1, int a2)` — genotype index of the unordered allele pair.
pub fn genotype(a1: i32, a2: i32) -> i32 {
    if a1 <= a2 {
        assert!(a1 >= 0, "allele < 0: {a1} {a2}");
        (a2 * (a2 + 1)) / 2 + a1
    } else {
        assert!(a2 >= 0, "allele < 0: {a1} {a2}");
        (a1 * (a1 + 1)) / 2 + a2
    }
}

fn check_samples(samples: &Samples, recs: &[Rc<dyn GTRec>]) {
    for rec in recs {
        assert!(rec.samples() == samples, "inconsistent samples");
    }
}

fn check_markers_and_samples(markers: &Markers, samples: &Samples, recs: &[Rc<dyn GTRec>]) {
    for (j, rec) in recs.iter().enumerate() {
        assert!(
            rec.marker() == markers.marker(j as i32),
            "inconsistent markers"
        );
        assert!(rec.samples() == samples, "inconsistent samples");
    }
}

fn markers_of(recs: &[Rc<dyn GTRec>]) -> Markers {
    let markers: Vec<Marker> = recs.iter().map(|r| r.marker().clone()).collect();
    Markers::create(markers)
}

fn is_ref_data(recs: &[Rc<dyn GTRec>]) -> bool {
    recs.iter().all(|r| r.is_phased())
}

impl BasicGT {
    /// `new BasicGT(GTRec[] recs)`.
    pub fn new(recs: Vec<Rc<dyn GTRec>>) -> Self {
        let samples = recs[0].samples().clone();
        BasicGT::from_samples(samples, recs)
    }

    /// `new BasicGT(Samples samples, GTRec[] recs)`.
    pub fn from_samples(samples: Samples, recs: Vec<Rc<dyn GTRec>>) -> Self {
        check_samples(&samples, &recs);
        let markers = markers_of(&recs);
        let is_ref_data = is_ref_data(&recs);
        BasicGT {
            samples,
            markers,
            recs,
            is_ref_data,
        }
    }

    /// `new BasicGT(Markers markers, Samples samples, GTRec[] recs)`.
    pub fn from_markers_samples(
        markers: Markers,
        samples: Samples,
        recs: Vec<Rc<dyn GTRec>>,
    ) -> Self {
        check_markers_and_samples(&markers, &samples, &recs);
        let is_ref_data = is_ref_data(&recs);
        BasicGT {
            samples,
            markers,
            recs,
            is_ref_data,
        }
    }

    /// `BasicGT.restrict(BasicGT gt, int[] indices)` — restrict to the given marker indices
    /// (markers recomputed from `gt`'s markers).
    pub fn restrict_indices(gt: &BasicGT, indices: &[i32]) -> BasicGT {
        let restrict_markers = gt.markers.restrict_indices(indices);
        let restricted_recs: Vec<Rc<dyn GTRec>> = indices
            .iter()
            .map(|&j| gt.recs[j as usize].clone())
            .collect();
        BasicGT::from_markers_samples(restrict_markers, gt.samples.clone(), restricted_recs)
    }
}

impl GT for BasicGT {
    fn is_reversed(&self) -> bool {
        false
    }

    fn n_markers(&self) -> i32 {
        self.recs.len() as i32
    }

    fn marker(&self, marker_index: i32) -> &Marker {
        self.markers.marker(marker_index)
    }

    fn markers(&self) -> &Markers {
        &self.markers
    }

    fn n_haps(&self) -> i32 {
        2 * self.samples.size()
    }

    fn n_samples(&self) -> i32 {
        self.samples.size()
    }

    fn samples(&self) -> &Samples {
        &self.samples
    }

    fn is_phased(&self) -> bool {
        self.is_ref_data
    }

    fn allele(&self, marker: i32, hap: i32) -> i32 {
        self.recs[marker as usize].get(hap)
    }

    fn restrict(&self, restricted_markers: &Markers, indices: &[i32]) -> Box<dyn GT> {
        let restricted_recs: Vec<Rc<dyn GTRec>> = indices
            .iter()
            .map(|&j| self.recs[j as usize].clone())
            .collect();
        Box::new(BasicGT::from_markers_samples(
            restricted_markers.clone(),
            self.samples.clone(),
            restricted_recs,
        ))
    }

    fn restrict_range(&self, start: i32, end: i32) -> Box<dyn GT> {
        let restrict_markers = self.markers.restrict(start, end);
        let restrict_recs: Vec<Rc<dyn GTRec>> = (start..end)
            .map(|j| self.recs[j as usize].clone())
            .collect();
        Box::new(BasicGT::from_markers_samples(
            restrict_markers,
            self.samples.clone(),
            restrict_recs,
        ))
    }
}

impl fmt::Display for BasicGT {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[BasicGT: nMarkers={} nSamples={}",
            self.n_markers(),
            self.n_samples()
        )?;
        for rec in &self.recs {
            f.write_str(consts::NL)?;
            f.write_str(&to_vcf_rec(rec.as_ref()))?;
        }
        f.write_str("]")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{BasicGTRec, MarkerParser, VcfHeader, HEADER_PREFIX};

    fn header(n_dip: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n_dip {
            hdr.push_str(&format!("\tS{s}"));
        }
        let lines = vec!["##fileformat=VCFv4.2".to_string(), hdr];
        let is_dip = vec![true; n_dip];
        VcfHeader::new_accept_all("src", &lines, &is_dip)
    }

    fn rec(h: &VcfHeader, line: &str) -> Rc<dyn GTRec> {
        let p =
            super::super::VcfRecGTParser::new(h, line, &MarkerParser::new(true, true, true, true));
        Rc::new(BasicGTRec::from_parser(&p))
    }

    #[test]
    fn genotype_indices() {
        // Standard VCF genotype-index ordering.
        assert_eq!(genotype(0, 0), 0);
        assert_eq!(genotype(0, 1), 1);
        assert_eq!(genotype(1, 1), 2);
        assert_eq!(genotype(0, 2), 3);
        assert_eq!(genotype(1, 2), 4);
        assert_eq!(genotype(2, 2), 5);
        // Unordered: genotype(a1,a2)==genotype(a2,a1)
        assert_eq!(genotype(2, 1), genotype(1, 2));
    }

    #[test]
    fn basic_construction_and_allele_access() {
        let h = header(2);
        let recs = vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0"),
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|1\t0|0"),
        ];
        let gt = BasicGT::new(recs);
        assert_eq!(gt.n_markers(), 2);
        assert_eq!(gt.n_samples(), 2);
        assert_eq!(gt.n_haps(), 4);
        assert!(!gt.is_reversed());
        assert!(gt.is_phased()); // all phased -> ref data
        assert_eq!(gt.allele(0, 0), 0);
        assert_eq!(gt.allele(0, 1), 1);
        assert_eq!(gt.allele(1, 0), 1);
        assert_eq!(gt.allele(1, 3), 0);
    }

    #[test]
    fn restrict_range_and_indices() {
        let h = header(1);
        let recs = vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1"),
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|0"),
            rec(&h, "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|0"),
        ];
        let gt = BasicGT::new(recs);
        let sub = gt.restrict_range(1, 3);
        assert_eq!(sub.n_markers(), 2);
        assert_eq!(sub.marker(0).pos(), 200);
        assert_eq!(sub.allele(0, 0), 1);

        let sub2 = BasicGT::restrict_indices(&gt, &[0, 2]);
        assert_eq!(sub2.n_markers(), 2);
        assert_eq!(sub2.marker(0).pos(), 100);
        assert_eq!(sub2.marker(1).pos(), 300);
    }

    #[test]
    fn unphased_is_not_ref_data() {
        let h = header(1);
        let recs = vec![rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1")];
        let gt = BasicGT::new(recs);
        assert!(!gt.is_phased());
    }
}
