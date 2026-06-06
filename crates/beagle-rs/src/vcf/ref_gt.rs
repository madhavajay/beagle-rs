//! Port of `vcf/RefGT.java` — reference (phased, non-missing) genotypes for a list of
//! markers and samples, backed by one `RefGTRec` per marker. Records are shared via `Rc`
//! across `restrict` views.

use std::rc::Rc;

use super::{Marker, Markers, RefGTRec, Samples, GT};

/// Port of `vcf/RefGT.java`.
#[derive(Clone)]
pub struct RefGT {
    markers: Markers,
    samples: Samples,
    recs: Vec<Rc<dyn RefGTRec>>,
}

fn check_data_recs(recs: &[Rc<dyn RefGTRec>]) -> Samples {
    assert!(!recs.is_empty(), "Missing data in VCF file");
    let samples = recs[0].samples().clone();
    for (j, rec) in recs.iter().enumerate() {
        assert!(
            rec.samples() == &samples,
            "sample inconsistency at index {j}"
        );
        assert!(rec.is_phased(), "non-reference data at marker index {j}");
    }
    samples
}

fn check_data(markers: &Markers, samples: &Samples, recs: &[Rc<dyn RefGTRec>]) {
    assert!(
        markers.size() == recs.len() as i32,
        "markers.nMarkers()={} refVcfRecs.length={}",
        markers.size(),
        recs.len()
    );
    for (j, rec) in recs.iter().enumerate() {
        assert!(
            rec.samples() == samples,
            "sample inconsistency at index {j}"
        );
        assert!(
            rec.marker() == markers.marker(j as i32),
            "marker inconsistency at index {j}"
        );
        assert!(rec.is_phased(), "non-reference data at marker index {j}");
    }
}

impl RefGT {
    /// `new RefGT(Markers, Samples, RefGTRec[])`.
    pub fn new(markers: Markers, samples: Samples, ref_vcf_recs: Vec<Rc<dyn RefGTRec>>) -> Self {
        check_data(&markers, &samples, &ref_vcf_recs);
        RefGT {
            markers,
            samples,
            recs: ref_vcf_recs,
        }
    }

    /// `new RefGT(RefGTRec[])`.
    pub fn from_recs(ref_vcf_recs: Vec<Rc<dyn RefGTRec>>) -> Self {
        let samples = check_data_recs(&ref_vcf_recs);
        let ma: Vec<Marker> = ref_vcf_recs.iter().map(|r| r.marker().clone()).collect();
        let markers = Markers::create(ma);
        RefGT {
            markers,
            samples,
            recs: ref_vcf_recs,
        }
    }

    /// `RefGT.restrict(RefGT refGT, int[] indices)` — select the given strictly increasing
    /// marker indices. Note (quirk vcf-5): Java passes `refGT.markers` (the FULL marker
    /// list) rather than a restricted list, so the resulting `RefGT` only passes its own
    /// `checkData` when `indices.length == refGT.nMarkers()`. This static method is dead in
    /// Beagle (only the instance `restrict` overloads are called); preserved faithfully.
    pub fn restrict_indices(ref_gt: &RefGT, indices: &[i32]) -> RefGT {
        let rra = select_strictly_increasing(&ref_gt.recs, indices);
        RefGT::new(ref_gt.markers.clone(), ref_gt.samples.clone(), rra)
    }

    /// `get(int marker)` — the `RefGTRec` for the marker.
    pub fn get(&self, marker: i32) -> Rc<dyn RefGTRec> {
        self.recs[marker as usize].clone()
    }

    /// `restrict(Markers, int[])` returning the concrete `RefGT` (the inherent form of the
    /// `GT::restrict` override).
    pub fn restrict_to_ref(&self, markers: &Markers, indices: &[i32]) -> RefGT {
        let rra = select_strictly_increasing(&self.recs, indices);
        RefGT::new(markers.clone(), self.samples.clone(), rra)
    }
}

fn select_strictly_increasing(recs: &[Rc<dyn RefGTRec>], indices: &[i32]) -> Vec<Rc<dyn RefGTRec>> {
    let mut rra = Vec::with_capacity(indices.len());
    for (j, &idx) in indices.iter().enumerate() {
        assert!(!(j > 0 && idx <= indices[j - 1]), "{idx}");
        rra.push(recs[idx as usize].clone());
    }
    rra
}

impl GT for RefGT {
    fn is_reversed(&self) -> bool {
        false
    }

    fn n_markers(&self) -> i32 {
        self.markers.size()
    }

    fn marker(&self, marker: i32) -> &Marker {
        self.markers.marker(marker)
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
        true
    }

    fn allele(&self, marker: i32, haplotype: i32) -> i32 {
        self.recs[marker as usize].get(haplotype)
    }

    fn restrict(self: Rc<Self>, markers: &Markers, indices: &[i32]) -> Rc<dyn GT> {
        Rc::new(self.restrict_to_ref(markers, indices))
    }

    fn restrict_range(self: Rc<Self>, start: i32, end: i32) -> Rc<dyn GT> {
        let restrict_markers = self.markers.restrict(start, end);
        let restrict_recs: Vec<Rc<dyn RefGTRec>> = (start..end)
            .map(|j| self.recs[j as usize].clone())
            .collect();
        Rc::new(RefGT::new(
            restrict_markers,
            self.samples.clone(),
            restrict_recs,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{allele_ref_gt_rec_from_components, MarkerParser};

    fn rec(
        line: &str,
        samples: &Samples,
        allele_to_haps: Vec<Option<Vec<i32>>>,
    ) -> Rc<dyn RefGTRec> {
        let marker = Marker::instance(line, &MarkerParser::new(true, true, true, true));
        Rc::from(allele_ref_gt_rec_from_components(
            marker,
            samples.clone(),
            allele_to_haps,
        ))
    }

    #[test]
    fn construction_and_allele_access() {
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        // marker0: 4 haps, allele 1 carried by haps [1,2] (major allele 0 -> None)
        let r0 = rec(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &samples,
            vec![None, Some(vec![1, 2])],
        );
        // marker1: allele 1 carried by hap [3]
        let r1 = rec(
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0",
            &samples,
            vec![None, Some(vec![3])],
        );
        let gt = RefGT::from_recs(vec![r0, r1]);
        assert_eq!(gt.n_markers(), 2);
        assert_eq!(gt.n_haps(), 4);
        assert!(gt.is_phased());
        // marker0 alleles per hap: [0,1,1,0]
        assert_eq!(
            (0..4).map(|h| gt.allele(0, h)).collect::<Vec<_>>(),
            vec![0, 1, 1, 0]
        );
        // marker1 alleles per hap: [0,0,0,1]
        assert_eq!(
            (0..4).map(|h| gt.allele(1, h)).collect::<Vec<_>>(),
            vec![0, 0, 0, 1]
        );
    }

    #[test]
    fn restrict_views() {
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        let mk = |p: i32| {
            rec(
                &format!("chr1\t{p}\t.\tA\tC\t.\tPASS\t.\tGT\t0|1"),
                &samples,
                vec![None, Some(vec![1])],
            )
        };
        let gt = Rc::new(RefGT::from_recs(vec![mk(100), mk(200), mk(300)]));
        let sub = gt.clone().restrict_range(1, 3);
        assert_eq!(sub.n_markers(), 2);
        assert_eq!(sub.marker(0).pos(), 200);

        // Instance restrict(markers, indices): caller supplies the restricted markers.
        let restricted_markers = gt.markers.restrict_indices(&[0, 2]);
        let sub2 = gt.clone().restrict(&restricted_markers, &[0, 2]);
        assert_eq!(sub2.n_markers(), 2);
        assert_eq!(sub2.marker(1).pos(), 300);
    }

    #[test]
    fn static_restrict_preserves_full_markers_quirk() {
        // vcf-5: static restrict passes the FULL markers, so it only validates when
        // indices covers all markers (here an identity permutation of the 2 markers).
        let samples = Samples::new(&["s0".into()], &[true]);
        let mk = |p: i32| {
            rec(
                &format!("chr1\t{p}\t.\tA\tC\t.\tPASS\t.\tGT\t0|1"),
                &samples,
                vec![None, Some(vec![1])],
            )
        };
        let gt = RefGT::from_recs(vec![mk(100), mk(200)]);
        let same = RefGT::restrict_indices(&gt, &[0, 1]);
        assert_eq!(same.n_markers(), 2);
        assert_eq!(same.marker(0).pos(), 100);
    }

    #[test]
    #[should_panic(expected = "0")]
    fn restrict_rejects_non_increasing_indices() {
        let samples = Samples::new(&["s0".into()], &[true]);
        let r = rec(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &samples,
            vec![None, Some(vec![1])],
        );
        let gt = RefGT::from_recs(vec![r]);
        // duplicate/non-increasing -> panic
        let _ = RefGT::restrict_indices(&gt, &[0, 0]);
    }
}
