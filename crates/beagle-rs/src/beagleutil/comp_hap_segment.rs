//! Port of `beagleutil/CompHapSegment.java` — a copied haplotype segment in a
//! composite reference haplotype. Ordered by `lastIbsStep` (matching `compareTo`).

use std::cmp::Ordering;

/// Port of `beagleutil/CompHapSegment.java`.
#[derive(Clone)]
pub struct CompHapSegment {
    hap: i32,
    start_marker: i32,
    last_ibs_step: i32,
    comp_hap_index: i32,
}

impl CompHapSegment {
    /// `new CompHapSegment(hap, startMarker, ibsStep, compHapIndex)`.
    pub fn new(hap: i32, start_marker: i32, ibs_step: i32, comp_hap_index: i32) -> Self {
        CompHapSegment {
            hap,
            start_marker,
            last_ibs_step: ibs_step,
            comp_hap_index,
        }
    }

    /// `updateSegment(hap, startMarker, lastIbsStep)`.
    pub fn update_segment(&mut self, hap: i32, start_marker: i32, last_ibs_step: i32) {
        self.hap = hap;
        self.start_marker = start_marker;
        self.last_ibs_step = last_ibs_step;
    }

    /// `setLastIbsStep(ibsStep)`.
    pub fn set_last_ibs_step(&mut self, ibs_step: i32) {
        self.last_ibs_step = ibs_step;
    }

    /// `hap()`.
    pub fn hap(&self) -> i32 {
        self.hap
    }

    /// `startMarker()`.
    pub fn start_marker(&self) -> i32 {
        self.start_marker
    }

    /// `lastIbsStep()`.
    pub fn last_ibs_step(&self) -> i32 {
        self.last_ibs_step
    }

    /// `compHapIndex()`.
    pub fn comp_hap_index(&self) -> i32 {
        self.comp_hap_index
    }
}

// `compareTo` orders by lastIbsStep only; Eq/Ord here mirror that exactly.
impl PartialEq for CompHapSegment {
    fn eq(&self, other: &Self) -> bool {
        self.last_ibs_step == other.last_ibs_step
    }
}
impl Eq for CompHapSegment {}
impl Ord for CompHapSegment {
    fn cmp(&self, other: &Self) -> Ordering {
        self.last_ibs_step.cmp(&other.last_ibs_step)
    }
}
impl PartialOrd for CompHapSegment {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordered_by_last_ibs_step() {
        let a = CompHapSegment::new(1, 10, 5, 100);
        let mut b = CompHapSegment::new(2, 20, 8, 200);
        assert!(a < b);
        b.set_last_ibs_step(5);
        assert_eq!(a.cmp(&b), Ordering::Equal);
        b.update_segment(9, 90, 3);
        assert_eq!(b.hap(), 9);
        assert_eq!(b.start_marker(), 90);
        assert!(b < a);
        assert_eq!(a.comp_hap_index(), 100);
    }
}
