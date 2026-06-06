//! Port of `vcf/MarkerMap.java` — genetic map positions and inter-marker genetic
//! distances for a sequence of loci.
//!
//! Parity note: `pRecomb` uses `expm1`. Java calls `Math.expm1`; the Rust port uses
//! `f64::exp_m1`. Both compute `e^x - 1` in IEEE-754 double precision but the underlying
//! libm implementations may differ by up to ~1 ULP. This is the one transcendental in this
//! file and is flagged for the end-to-end parity harness (see docs/known-quirks.md).

use crate::blbutil::{DoubleArray, FloatArray};
use crate::ints::IntArray;

use super::{gen_pos_markers_min_dist, GeneticMap, Markers};

/// Port of `vcf/MarkerMap.java`.
pub struct MarkerMap {
    gen_pos: DoubleArray,
    gen_dist: FloatArray,
}

fn gen_dist(gen_pos: &[f64]) -> FloatArray {
    let mut da = vec![0.0f32; gen_pos.len()];
    for j in 1..da.len() {
        da[j] = (gen_pos[j] - gen_pos[j - 1]) as f32;
        assert!(
            da[j] > 0.0,
            "Nonpositive genetic distance: dist[{j}]={}",
            da[j]
        );
    }
    FloatArray::from_floats(&da)
}

impl MarkerMap {
    fn from_gen_pos(g_pos: Vec<f64>) -> Self {
        let gen_dist = gen_dist(&g_pos);
        MarkerMap {
            gen_pos: DoubleArray::from_doubles(&g_pos),
            gen_dist,
        }
    }

    /// `MarkerMap.create(GeneticMap genMap, Markers markers)`.
    pub fn create(gen_map: &dyn GeneticMap, markers: &Markers) -> MarkerMap {
        let mean_gen_diff = mean_single_base_gen_dist(gen_map, markers);
        MarkerMap::from_gen_pos(gen_pos_markers_min_dist(gen_map, mean_gen_diff, markers))
    }

    /// `MarkerMap.create(GeneticMap genMap, double minGenDist, Markers markers)`.
    pub fn create_min_dist(
        gen_map: &dyn GeneticMap,
        min_gen_dist: f64,
        markers: &Markers,
    ) -> MarkerMap {
        MarkerMap::from_gen_pos(gen_pos_markers_min_dist(gen_map, min_gen_dist, markers))
    }

    /// `restrict(int[] indices)` — distinct, strictly increasing marker indices.
    pub fn restrict(&self, indices: &[i32]) -> MarkerMap {
        let mut g_pos = vec![0.0f64; indices.len()];
        g_pos[0] = self.gen_pos.get(indices[0]);
        for j in 1..indices.len() {
            assert!(indices[j] > indices[j - 1], "{}", indices[j]);
            g_pos[j] = self.gen_pos.get(indices[j]);
        }
        MarkerMap::from_gen_pos(g_pos)
    }

    /// `restrict(IntArray indices)`.
    pub fn restrict_int_array(&self, indices: &dyn IntArray) -> MarkerMap {
        let n = indices.size() as usize;
        let mut g_pos = vec![0.0f64; n];
        g_pos[0] = self.gen_pos.get(indices.get(0));
        #[allow(clippy::needless_range_loop)] // mirrors Java indexed loop reading j and j-1
        for j in 1..n {
            assert!(
                indices.get(j as i32) > indices.get(j as i32 - 1),
                "{}",
                indices.get(j as i32)
            );
            g_pos[j] = self.gen_pos.get(indices.get(j as i32));
        }
        MarkerMap::from_gen_pos(g_pos)
    }

    /// `genPos()` — genetic map position of each marker.
    pub fn gen_pos(&self) -> &DoubleArray {
        &self.gen_pos
    }

    /// `genDist()` — genetic distance to the previous marker (`0.0` at index 0).
    pub fn gen_dist(&self) -> &FloatArray {
        &self.gen_dist
    }

    /// `pRecomb(float recombIntensity)` — recombination probability per inter-marker gap.
    pub fn p_recomb(&self, recomb_intensity: f32) -> FloatArray {
        assert!(
            recomb_intensity > 0.0 && recomb_intensity.is_finite(),
            "{recomb_intensity}"
        );
        let c = -(recomb_intensity as f64);
        let p: Vec<f64> = (0..self.gen_dist.size())
            .map(|m| -((c * self.gen_dist.get(m) as f64).exp_m1()))
            .collect();
        FloatArray::from_doubles(&p)
    }
}

/// `MarkerMap.meanSingleBaseGenDist(GeneticMap genMap, Markers markers)`.
pub fn mean_single_base_gen_dist(gen_map: &dyn GeneticMap, markers: &Markers) -> f64 {
    let a = markers.marker(0);
    let b = markers.marker(markers.size() - 1);
    assert!(a.chrom_index() == b.chrom_index(), "inconsistent data");
    assert!(
        a.pos() != b.pos(),
        "Window has only one position: CHROM={} POS={}",
        a.chrom(),
        a.pos()
    );
    let mean_single_base_dist = (gen_map.gen_pos_marker(b) - gen_map.gen_pos_marker(a)).abs()
        / ((b.pos() - a.pos()).abs() as f64);
    // require >= 0.01 * mean human single-base genetic distance
    mean_single_base_dist.max(1e-8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{Marker, MarkerParser, PositionMap};

    fn markers(positions: &[i32]) -> Markers {
        let mp = MarkerParser::new(true, true, true, true);
        let ms: Vec<Marker> = positions
            .iter()
            .map(|p| Marker::instance(&format!("chr1\t{p}\t.\tA\tC\t.\tPASS\t.\tGT\t0|0"), &mp))
            .collect();
        Markers::create(ms)
    }

    #[test]
    fn gen_pos_and_dist_with_position_map() {
        // PositionMap(scale): genPos = pos * scale. scale 1e-6 -> 1bp = 1e-6 cM.
        let gm = PositionMap::new(1e-6);
        let ms = markers(&[1_000_000, 2_000_000, 3_500_000]);
        let mm = MarkerMap::create_min_dist(&gm, 0.0, &ms);
        // genPos
        assert!((mm.gen_pos().get(0) - 1.0).abs() < 1e-12);
        assert!((mm.gen_pos().get(1) - 2.0).abs() < 1e-12);
        assert!((mm.gen_pos().get(2) - 3.5).abs() < 1e-12);
        // genDist: [0, 1.0, 1.5]
        assert_eq!(mm.gen_dist().get(0), 0.0);
        assert!((mm.gen_dist().get(1) - 1.0).abs() < 1e-6);
        assert!((mm.gen_dist().get(2) - 1.5).abs() < 1e-6);
    }

    #[test]
    fn p_recomb_is_one_minus_exp() {
        let gm = PositionMap::new(1e-6);
        let ms = markers(&[1_000_000, 2_000_000]);
        let mm = MarkerMap::create_min_dist(&gm, 0.0, &ms);
        // gap dist = 1.0 cM, intensity = 0.5 -> p = -expm1(-0.5*1.0) = 1 - e^-0.5
        let p = mm.p_recomb(0.5);
        let expected = 1.0f64 - (-0.5f64).exp();
        assert!((p.get(1) as f64 - expected).abs() < 1e-6);
        assert_eq!(p.get(0), 0.0); // genDist[0]=0 -> -expm1(0)=0
    }

    #[test]
    fn restrict_subsets_positions() {
        let gm = PositionMap::new(1e-6);
        let ms = markers(&[1_000_000, 2_000_000, 3_000_000, 4_000_000]);
        let mm = MarkerMap::create_min_dist(&gm, 0.0, &ms);
        let sub = mm.restrict(&[0, 2, 3]);
        assert!((sub.gen_pos().get(0) - 1.0).abs() < 1e-12);
        assert!((sub.gen_pos().get(1) - 3.0).abs() < 1e-12);
        assert!((sub.gen_pos().get(2) - 4.0).abs() < 1e-12);
        // dist between restricted positions 1.0->3.0 = 2.0
        assert!((sub.gen_dist().get(1) - 2.0).abs() < 1e-6);
    }
}
