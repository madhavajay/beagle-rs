//! Port of `vcf/GeneticMap.java` — genetic-map interface + static helpers. The
//! `geneticMap(file, chromInt)` factory is completed once `PlinkGenMap` is ported; the
//! `PositionMap` default (coordinate × scale) is available now.

use super::{Marker, Markers};

/// Port of the `vcf/GeneticMap.java` interface.
pub trait GeneticMap {
    /// `basePos(int chrom, double geneticPosition)`.
    fn base_pos(&self, chrom: i32, genetic_position: f64) -> i32;

    /// `genPos(Marker marker)`.
    fn gen_pos_marker(&self, marker: &Marker) -> f64;

    /// `genPos(int chrom, int basePosition)`.
    fn gen_pos(&self, chrom: i32, base_position: i32) -> f64;
}

/// `GeneticMap.genPos(genMap, markers)` — map position of each marker.
pub fn gen_pos_markers(gen_map: &dyn GeneticMap, markers: &Markers) -> Vec<f64> {
    assert!(
        markers.marker(0).chrom_index() == markers.marker(markers.size() - 1).chrom_index(),
        "inconsistent data"
    );
    (0..markers.size())
        .map(|j| gen_map.gen_pos_marker(markers.marker(j)))
        .collect()
}

/// `GeneticMap.genPos(genMap, minGenDist, markers)` — map positions enforcing a minimum
/// cM distance between successive markers.
pub fn gen_pos_markers_min_dist(
    gen_map: &dyn GeneticMap,
    min_gen_dist: f64,
    markers: &Markers,
) -> Vec<f64> {
    assert!(
        markers.marker(0).chrom_index() == markers.marker(markers.size() - 1).chrom_index(),
        "inconsistent data"
    );
    assert!(min_gen_dist.is_finite(), "{}", min_gen_dist);
    let n = markers.size() as usize;
    let mut gen_pos = vec![0.0f64; n];
    gen_pos[0] = gen_map.gen_pos_marker(markers.marker(0));
    let mut last_map_pos = gen_pos[0];
    for j in 1..n {
        let map_pos = gen_map.gen_pos_marker(markers.marker(j as i32));
        let dist = (map_pos - last_map_pos).max(min_gen_dist);
        gen_pos[j] = gen_pos[j - 1] + dist;
        last_map_pos = map_pos;
    }
    gen_pos
}
