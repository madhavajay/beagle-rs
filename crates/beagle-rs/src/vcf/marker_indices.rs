//! Port of `vcf/MarkerIndices.java` — overlap/splice indices for a marker window plus the
//! bidirectional mapping between all-marker indices and target-marker indices.

use crate::ints::java_binary_search;

/// Port of `vcf/MarkerIndices.java`.
#[derive(Clone)]
pub struct MarkerIndices {
    prev_splice: i32,
    overlap_end: i32,
    overlap_start: i32,
    next_splice: i32,
    targ_marker_to_marker: Vec<i32>,
    marker_to_targ_marker: Vec<i32>,
    targ_prev_splice: i32,
    targ_overlap_end: i32,
    targ_overlap_start: i32,
    targ_next_splice: i32,
}

/// First target marker on or after `marker` (insertion point in the sorted map).
fn targ_index(targ_marker_to_marker: &[i32], marker: i32) -> i32 {
    let ins_pt = java_binary_search(
        targ_marker_to_marker,
        0,
        targ_marker_to_marker.len() as i32,
        marker,
    );
    if ins_pt < 0 {
        -ins_pt - 1
    } else {
        ins_pt
    }
}

fn targ_marker_to_marker_from(in_targ: &[bool]) -> Vec<i32> {
    // Java sizes an IntList with `1 + inTarg.length>>6` == `(1 + len) >> 6` (capacity only).
    let mut il = Vec::with_capacity((1 + in_targ.len()) >> 6);
    for (j, &t) in in_targ.iter().enumerate() {
        if t {
            il.push(j as i32);
        }
    }
    il
}

fn marker_to_targ_marker_from(targ_marker_to_marker: &[i32], n_markers: usize) -> Vec<i32> {
    let mut ia = vec![-1i32; n_markers];
    for (j, &m) in targ_marker_to_marker.iter().enumerate() {
        ia[m as usize] = j as i32;
    }
    ia
}

impl MarkerIndices {
    /// `new MarkerIndices(boolean[] inTarg, int overlapEnd, int overlapStart)`.
    pub fn from_in_targ(in_targ: &[bool], overlap_end: i32, overlap_start: i32) -> Self {
        let n = in_targ.len() as i32;
        assert!(overlap_end >= 0 && overlap_end <= n, "{overlap_end}");
        assert!(overlap_start >= 0 && overlap_start <= n, "{overlap_start}");
        let prev_splice = overlap_end >> 1;
        let next_splice = (n + overlap_start) >> 1; // values non-negative: >>> == >>
        let targ_marker_to_marker = targ_marker_to_marker_from(in_targ);
        let marker_to_targ_marker =
            marker_to_targ_marker_from(&targ_marker_to_marker, in_targ.len());
        let targ_prev_splice = targ_index(&targ_marker_to_marker, prev_splice);
        let targ_overlap_end = targ_index(&targ_marker_to_marker, overlap_end);
        let targ_overlap_start = targ_index(&targ_marker_to_marker, overlap_start);
        let targ_next_splice = targ_index(&targ_marker_to_marker, next_splice);
        MarkerIndices {
            prev_splice,
            overlap_end,
            overlap_start,
            next_splice,
            targ_marker_to_marker,
            marker_to_targ_marker,
            targ_prev_splice,
            targ_overlap_end,
            targ_overlap_start,
            targ_next_splice,
        }
    }

    /// `new MarkerIndices(int overlapEnd, int overlapStart, int nMarkers)` — identity map
    /// (every marker is a target marker).
    pub fn from_counts(overlap_end: i32, overlap_start: i32, n_markers: i32) -> Self {
        assert!(n_markers >= 0, "{n_markers}");
        assert!(
            overlap_end >= 0 && overlap_end <= n_markers,
            "{overlap_end}"
        );
        assert!(
            overlap_start >= 0 && overlap_start <= n_markers,
            "{overlap_start}"
        );
        let prev_splice = overlap_end >> 1;
        let next_splice = (n_markers + overlap_start) >> 1;
        let identity: Vec<i32> = (0..n_markers).collect();
        MarkerIndices {
            prev_splice,
            overlap_end,
            overlap_start,
            next_splice,
            targ_marker_to_marker: identity.clone(),
            marker_to_targ_marker: identity,
            targ_prev_splice: prev_splice,
            targ_overlap_end: overlap_end,
            targ_overlap_start: overlap_start,
            targ_next_splice: next_splice,
        }
    }

    /// `nMarkers()`.
    pub fn n_markers(&self) -> i32 {
        self.marker_to_targ_marker.len() as i32
    }

    /// `nTargMarkers()`.
    pub fn n_targ_markers(&self) -> i32 {
        self.targ_marker_to_marker.len() as i32
    }

    /// `prevSplice()`.
    pub fn prev_splice(&self) -> i32 {
        self.prev_splice
    }

    /// `overlapEnd()`.
    pub fn overlap_end(&self) -> i32 {
        self.overlap_end
    }

    /// `overlapStart()`.
    pub fn overlap_start(&self) -> i32 {
        self.overlap_start
    }

    /// `prevTargSplice()`.
    pub fn prev_targ_splice(&self) -> i32 {
        self.targ_prev_splice
    }

    /// `nextSplice()`.
    pub fn next_splice(&self) -> i32 {
        self.next_splice
    }

    /// `targOverlapEnd()`.
    pub fn targ_overlap_end(&self) -> i32 {
        self.targ_overlap_end
    }

    /// `targOverlapStart()`.
    pub fn targ_overlap_start(&self) -> i32 {
        self.targ_overlap_start
    }

    /// `nextTargSplice()`.
    pub fn next_targ_splice(&self) -> i32 {
        self.targ_next_splice
    }

    /// `targMarkerToMarker(int targetMarker)`.
    pub fn targ_marker_to_marker_at(&self, target_marker: i32) -> i32 {
        self.targ_marker_to_marker[target_marker as usize]
    }

    /// `targMarkerToMarker()` — the full map (length `nTargMarkers`).
    pub fn targ_marker_to_marker(&self) -> Vec<i32> {
        self.targ_marker_to_marker.clone()
    }

    /// `markerToTargMarker(int marker)` — target index, or -1 if not in the target data.
    pub fn marker_to_targ_marker_at(&self, marker: i32) -> i32 {
        self.marker_to_targ_marker[marker as usize]
    }

    /// `markerToTargMarker()` — the full map (length `nMarkers`).
    pub fn marker_to_targ_marker(&self) -> Vec<i32> {
        self.marker_to_targ_marker.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn in_targ_mapping() {
        // 5 markers; targets at indices 0, 2, 4
        let in_targ = [true, false, true, false, true];
        let mi = MarkerIndices::from_in_targ(&in_targ, 2, 4);
        assert_eq!(mi.n_markers(), 5);
        assert_eq!(mi.n_targ_markers(), 3);
        assert_eq!(mi.targ_marker_to_marker(), vec![0, 2, 4]);
        assert_eq!(mi.marker_to_targ_marker(), vec![0, -1, 1, -1, 2]);
        assert_eq!(mi.targ_marker_to_marker_at(1), 2);
        assert_eq!(mi.marker_to_targ_marker_at(3), -1);
        // splices: prevSplice = 2>>1 = 1; nextSplice = (5+4)>>1 = 4
        assert_eq!(mi.prev_splice(), 1);
        assert_eq!(mi.next_splice(), 4);
        assert_eq!(mi.overlap_end(), 2);
        assert_eq!(mi.overlap_start(), 4);
        // targIndex(prevSplice=1) -> first target >= 1 is index 1 (marker 2)
        assert_eq!(mi.prev_targ_splice(), 1);
        // targIndex(overlapEnd=2) -> target index 1
        assert_eq!(mi.targ_overlap_end(), 1);
        // targIndex(overlapStart=4) -> target index 2
        assert_eq!(mi.targ_overlap_start(), 2);
        // targIndex(nextSplice=4) -> target index 2
        assert_eq!(mi.next_targ_splice(), 2);
    }

    #[test]
    fn identity_map_constructor() {
        let mi = MarkerIndices::from_counts(2, 6, 8);
        assert_eq!(mi.n_markers(), 8);
        assert_eq!(mi.n_targ_markers(), 8);
        assert_eq!(mi.targ_marker_to_marker(), (0..8).collect::<Vec<_>>());
        assert_eq!(mi.marker_to_targ_marker(), (0..8).collect::<Vec<_>>());
        assert_eq!(mi.prev_splice(), 1); // 2>>1
        assert_eq!(mi.next_splice(), 7); // (8+6)>>1
        assert_eq!(mi.prev_targ_splice(), 1);
        assert_eq!(mi.next_targ_splice(), 7);
        assert_eq!(mi.targ_overlap_end(), 2);
        assert_eq!(mi.targ_overlap_start(), 6);
    }
}
