//! Port of `vcf/Steps.java` — a partition of a marker list into consecutive "steps", each
//! spanning at least `minStep` cM (by first-marker genetic position).

use std::rc::Rc;

use super::MarkerMap;

/// Port of `vcf/Steps.java`.
pub struct Steps {
    map: Rc<MarkerMap>,
    step_ends: Vec<i32>,
}

fn step_ends(map: &MarkerMap, min_step: f64) -> Vec<i32> {
    let gen_pos = map.gen_pos();
    let n_markers = gen_pos.size();
    let mut indices: Vec<i32> = Vec::with_capacity((n_markers >> 1).max(0) as usize);
    let mut end = 0;
    while end < n_markers {
        let min_gen_pos = gen_pos.get(end) + min_step;
        end += 1;
        while end < gen_pos.size() && gen_pos.get(end) < min_gen_pos {
            end += 1;
        }
        indices.push(end);
    }
    indices
}

impl Steps {
    /// `new Steps(MarkerMap map, float minStep)`.
    pub fn new(map: Rc<MarkerMap>, min_step: f32) -> Self {
        assert!(min_step > 0.0 && min_step.is_finite(), "{min_step}");
        let step_ends = step_ends(&map, min_step as f64);
        Steps { map, step_ends }
    }

    /// `size()` — number of steps.
    pub fn size(&self) -> i32 {
        self.step_ends.len() as i32
    }

    /// `start(int step)` — first marker index of the step.
    pub fn start(&self, step: i32) -> i32 {
        if step == 0 {
            0
        } else {
            self.step_ends[(step - 1) as usize]
        }
    }

    /// `end(int step)` — last marker index (exclusive) of the step.
    pub fn end(&self, step: i32) -> i32 {
        self.step_ends[step as usize]
    }

    /// `map()` — the marker map.
    pub fn map(&self) -> &MarkerMap {
        &self.map
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{Marker, MarkerParser, Markers, PositionMap};

    fn markers(positions: &[i32]) -> Markers {
        let mp = MarkerParser::new(true, true, true, true);
        let ms: Vec<Marker> = positions
            .iter()
            .map(|p| Marker::instance(&format!("chr1\t{p}\t.\tA\tC\t.\tPASS\t.\tGT\t0|0"), &mp))
            .collect();
        Markers::create(ms)
    }

    fn map(positions: &[i32]) -> Rc<MarkerMap> {
        let gm = PositionMap::new(1e-6); // 1bp = 1e-6 cM
        Rc::new(MarkerMap::create_min_dist(&gm, 0.0, &markers(positions)))
    }

    #[test]
    fn partitions_by_min_step() {
        // positions in cM: 0.0, 0.5, 1.0, 1.5, 2.0 (scale 1e-6 over Mb-spaced bp)
        let mm = map(&[0, 500_000, 1_000_000, 1_500_000, 2_000_000]);
        // minStep 1.0 cM: step0 from genPos 0 needs >=1.0 -> ends at first marker with genPos>=1.0
        let steps = Steps::new(mm, 1.0);
        // step0: start 0; minGenPos=1.0; advance while <1.0: markers at 0.5 (<1.0) included,
        //   marker at 1.0 stops -> end=2. step1: from 1.0, minGenPos=2.0; 1.5<2.0 -> end=4;
        //   marker 2.0 stops -> end=4. step2: from 2.0, minGenPos=3.0; end=5 (past last).
        assert_eq!(steps.size(), 3);
        assert_eq!((steps.start(0), steps.end(0)), (0, 2));
        assert_eq!((steps.start(1), steps.end(1)), (2, 4));
        assert_eq!((steps.start(2), steps.end(2)), (4, 5));
    }

    #[test]
    fn small_step_makes_singletons() {
        let mm = map(&[0, 1_000_000, 2_000_000]);
        // tiny minStep -> each marker its own step
        let steps = Steps::new(mm, 0.001);
        assert_eq!(steps.size(), 3);
        assert_eq!((steps.start(0), steps.end(0)), (0, 1));
        assert_eq!((steps.start(2), steps.end(2)), (2, 3));
    }
}
