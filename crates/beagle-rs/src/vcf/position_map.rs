//! Port of `vcf/PositionMap.java` — a genetic map that multiplies base position by a
//! fixed scale factor (the default map when no PLINK map file is supplied; Beagle uses
//! `1 cM = 1 Mb`, i.e. scale `1e-6`).

use crate::beagleutil::ChromIds;

use super::{GeneticMap, Marker};

/// Port of `vcf/PositionMap.java`.
pub struct PositionMap {
    scale_factor: f64,
    inv_scale_factor: f64,
}

impl PositionMap {
    /// `new PositionMap(double scaleFactor)`.
    pub fn new(scale_factor: f64) -> Self {
        assert!(
            scale_factor.is_finite() && scale_factor > 0.0,
            "{}",
            scale_factor
        );
        PositionMap {
            scale_factor,
            inv_scale_factor: 1.0 / scale_factor,
        }
    }

    /// `scaleFactor()`.
    pub fn scale_factor(&self) -> f64 {
        self.scale_factor
    }
}

impl GeneticMap for PositionMap {
    fn base_pos(&self, chrom: i32, genetic_position: f64) -> i32 {
        assert!(
            chrom >= 0 && chrom < ChromIds::instance().size(),
            "{}",
            chrom
        );
        let pos = (genetic_position * self.inv_scale_factor).round() as i64;
        assert!(
            pos <= i32::MAX as i64,
            "An estimated base position exceeds the maximum integer value\n\
             Is the window parameter in cM units?"
        );
        pos as i32
    }

    fn gen_pos_marker(&self, marker: &Marker) -> f64 {
        self.scale_factor * marker.pos() as f64
    }

    fn gen_pos(&self, chrom: i32, base_position: i32) -> f64 {
        assert!(
            chrom >= 0 && chrom < ChromIds::instance().size(),
            "{}",
            chrom
        );
        self.scale_factor * base_position as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::beagleutil::ChromIds;
    use crate::vcf::MarkerParser;

    #[test]
    fn scale_map() {
        // ensure the chrom is registered so bounds checks pass
        let chrom = ChromIds::instance().get_index("chrPM");
        let pm = PositionMap::new(1e-6); // 1 cM = 1 Mb
        assert_eq!(pm.scale_factor(), 1e-6);
        assert!((pm.gen_pos(chrom, 2_000_000) - 2.0).abs() < 1e-12);
        assert_eq!(pm.base_pos(chrom, 2.0), 2_000_000);

        let marker = Marker::instance(
            "chrPM\t3000000\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &MarkerParser::new(true, true, true, true),
        );
        assert!((pm.gen_pos_marker(&marker) - 3.0).abs() < 1e-12);
    }

    #[test]
    #[should_panic]
    fn rejects_nonpositive_scale() {
        let _ = PositionMap::new(0.0);
    }
}
