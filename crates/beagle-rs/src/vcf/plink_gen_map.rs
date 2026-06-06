//! Port of `vcf/PlinkGenMap.java` — a genetic map from a PLINK map file (cM units),
//! with linear interpolation and a 5 cM end-extrapolation rule.

use crate::beagleutil::ChromIds;
use crate::blbutil::{binary_search_f64, Filter, InputIt, StringUtil};
use crate::ints::java_binary_search;
use std::path::Path;

use super::{GeneticMap, Marker};

const MIN_END_CM_DIST: f64 = 5.0;

/// Port of `vcf/PlinkGenMap.java`.
pub struct PlinkGenMap {
    base_pos: Vec<Vec<i32>>,
    gen_pos: Vec<Vec<f64>>,
}

fn fill_map_positions(list: &[String]) -> (Vec<i32>, Vec<f64>) {
    let n = list.len();
    let mut base_pos = vec![0i32; n];
    let mut gen_pos = vec![0.0f64; n];
    for j in 0..n {
        let fields = StringUtil::get_fields_ws(&list[j]);
        assert!(fields.len() == 4, "Map file format error: {}", list[j]);
        base_pos[j] = fields[3].parse().expect("parsable base position");
        gen_pos[j] = fields[2].parse().expect("parsable genetic position");
        assert!(
            gen_pos[j].is_finite(),
            "invalid map position: {}",
            gen_pos[j]
        );
        if j > 0 {
            assert!(
                base_pos[j] != base_pos[j - 1],
                "duplication position: {}",
                list[j]
            );
            assert!(
                !(base_pos[j] < base_pos[j - 1] || gen_pos[j] < gen_pos[j - 1]),
                "map positions not in ascending order: {}",
                list[j]
            );
        }
    }
    assert!(
        !(n > 0 && gen_pos[0] == gen_pos[n - 1]),
        "All loci in genetic map have the same genetic position [{}]: {}",
        gen_pos[0],
        list[0]
    );
    (base_pos, gen_pos)
}

impl PlinkGenMap {
    fn new(chrom_list: Vec<Vec<String>>) -> Self {
        let mut base_pos = Vec::with_capacity(chrom_list.len());
        let mut gen_pos = Vec::with_capacity(chrom_list.len());
        for list in &chrom_list {
            let (bp, gp) = fill_map_positions(list);
            base_pos.push(bp);
            gen_pos.push(gp);
        }
        PlinkGenMap { base_pos, gen_pos }
    }

    /// `PlinkGenMap.fromPlinkMapFile(File)`.
    pub fn from_plink_map_file(map_file: &Path) -> Self {
        PlinkGenMap::new(divide_by_chrom(map_file, &Filter::accept_all()))
    }

    /// `PlinkGenMap.fromPlinkMapFile(File, String chrom)`.
    pub fn from_plink_map_file_chrom(map_file: &Path, chrom: &str) -> Self {
        let filter = Filter::include([chrom.trim().to_string()]);
        PlinkGenMap::new(divide_by_chrom(map_file, &filter))
    }

    fn check_chrom_index(&self, chrom: i32) {
        assert!(
            chrom >= 0 && chrom < ChromIds::instance().size(),
            "{}",
            chrom
        );
        assert!(
            (chrom as usize) < self.base_pos.len() && !self.base_pos[chrom as usize].is_empty(),
            "missing genetic map for chromosome {}",
            ChromIds::instance().id(chrom)
        );
    }

    /// `nMapPositions(int chrom)`.
    pub fn n_map_positions(&self, chrom: i32) -> i32 {
        self.check_chrom_index(chrom);
        self.base_pos[chrom as usize].len() as i32
    }

    /// `index2BasePos(int chrom, int index)`.
    pub fn index2base_pos(&self, chrom: i32, index: i32) -> i32 {
        self.check_chrom_index(chrom);
        self.base_pos[chrom as usize][index as usize]
    }

    /// `index2GenPos(int chrom, int index)`.
    pub fn index2gen_pos(&self, chrom: i32, index: i32) -> f64 {
        self.check_chrom_index(chrom);
        self.gen_pos[chrom as usize][index as usize]
    }

    /// `closestIndex(int chrom, int basePosition)`. Note (vcf-3): the out-of-range
    /// branches compare against `basePos.length` (the chromosome count, outer array),
    /// not `basePos[chrom].length`, faithfully preserved from Java.
    pub fn closest_index(&self, chrom: i32, base_position: i32) -> i32 {
        self.check_chrom_index(chrom);
        let bp = &self.base_pos[chrom as usize];
        let map_index = java_binary_search(bp, 0, bp.len() as i32, base_position);
        if map_index >= 0 {
            return map_index;
        }
        let ins_pt = -map_index - 1;
        if ins_pt == 0 {
            0
        } else if ins_pt == self.base_pos.len() as i32 {
            // vcf-3: uses outer length (chromosome count), preserved
            self.base_pos.len() as i32 - 1
        } else {
            let dist_ins_pt = bp[ins_pt as usize] - base_position;
            let dist_ins_pt_m1 = base_position - bp[(ins_pt - 1) as usize];
            if dist_ins_pt <= dist_ins_pt_m1 {
                ins_pt
            } else {
                ins_pt - 1
            }
        }
    }
}

fn divide_by_chrom(map_file: &Path, chrom_filter: &Filter<String>) -> Vec<Vec<String>> {
    let mut chrom_list: Vec<Vec<String>> = Vec::new();
    for line in InputIt::from_gzip_file(map_file) {
        let fields = StringUtil::get_fields_ws_limit(&line, 4);
        if !fields.is_empty() {
            assert!(fields.len() >= 4, "Map file format error: {line}");
            let chrom = fields[0];
            if chrom_filter.accept(&chrom.to_string()) {
                let chrom_index = ChromIds::instance().get_index(chrom) as usize;
                while chrom_index >= chrom_list.len() {
                    chrom_list.push(Vec::new());
                }
                chrom_list[chrom_index].push(line.clone());
            }
        }
    }
    chrom_list
}

impl GeneticMap for PlinkGenMap {
    fn gen_pos_marker(&self, marker: &Marker) -> f64 {
        self.gen_pos(marker.chrom_index(), marker.pos())
    }

    fn gen_pos(&self, chrom: i32, base_position: i32) -> f64 {
        self.check_chrom_index(chrom);
        let bp = &self.base_pos[chrom as usize];
        let gp = &self.gen_pos[chrom as usize];
        let map_size_m1 = bp.len() as i32 - 1;
        let index = java_binary_search(bp, 0, bp.len() as i32, base_position);
        if index >= 0 {
            return gp[index as usize];
        }
        let ins_pt = -index - 1;
        let mut a_index = ins_pt - 1;
        let mut b_index = ins_pt;
        if a_index == map_size_m1 {
            let mut p = binary_search_f64(
                gp,
                0,
                gp.len() as i32,
                gp[map_size_m1 as usize] - MIN_END_CM_DIST,
            );
            if p < 0 {
                p = -p - 2;
            }
            a_index = p.max(0);
            b_index = map_size_m1;
        } else if b_index == 0 {
            let mut p = binary_search_f64(gp, 0, gp.len() as i32, gp[0] + MIN_END_CM_DIST);
            if p < 0 {
                p = -p - 1;
            }
            a_index = 0;
            b_index = p.min(map_size_m1);
        }
        let x = base_position;
        let a = bp[a_index as usize];
        let b = bp[b_index as usize];
        let fa = gp[a_index as usize];
        let fb = gp[b_index as usize];
        fa + (((x - a) as f64 / (b - a) as f64) * (fb - fa))
    }

    fn base_pos(&self, chrom: i32, genetic_position: f64) -> i32 {
        self.check_chrom_index(chrom);
        let bp = &self.base_pos[chrom as usize];
        let gp = &self.gen_pos[chrom as usize];
        let map_size_m1 = gp.len() as i32 - 1;
        let index = binary_search_f64(gp, 0, gp.len() as i32, genetic_position);
        if index >= 0 {
            return bp[index as usize];
        }
        let ins_pt = -index - 1;
        let mut a_index = ins_pt - 1;
        let mut b_index = ins_pt;
        if a_index == map_size_m1 {
            let mut p = binary_search_f64(
                gp,
                0,
                gp.len() as i32,
                gp[map_size_m1 as usize] - MIN_END_CM_DIST,
            );
            if p < 0 {
                p = -p - 2;
            }
            a_index = p.max(0);
            b_index = map_size_m1;
        } else if b_index == 0 {
            let mut p = binary_search_f64(gp, 0, gp.len() as i32, gp[0] + MIN_END_CM_DIST);
            if p < 0 {
                p = -p - 1;
            }
            a_index = 0;
            b_index = p.min(map_size_m1);
        }
        let x = genetic_position;
        let a = gp[a_index as usize];
        let b = gp[b_index as usize];
        let fa = bp[a_index as usize];
        let fb = bp[b_index as usize];
        let interp = fa as f64 + ((x - a) / (b - a)) * ((fb - fa) as f64);
        assert!(
            interp < i32::MAX as f64,
            "An estimated base position exceeds the maximum integer value\n\
             Are the window parameter and the genetic map in cM units?"
        );
        (interp + 0.5).floor() as i32 // Math.round
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_map(name: &str, contents: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let path = dir.join(name);
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    #[test]
    fn linear_interpolation() {
        // chrom id cM bp
        let map = "chrPLINK . 1.0 1000000\nchrPLINK . 2.0 2000000\nchrPLINK . 3.0 3000000\n";
        let path = write_map("beagle_rs_test_plink.map", map);
        let gm = PlinkGenMap::from_plink_map_file(&path);
        let chrom = ChromIds::instance().get_index("chrPLINK");
        assert_eq!(gm.n_map_positions(chrom), 3);
        // exact map points
        assert!((gm.gen_pos(chrom, 2_000_000) - 2.0).abs() < 1e-12);
        // interpolated midpoint
        assert!((gm.gen_pos(chrom, 1_500_000) - 1.5).abs() < 1e-12);
        // base position from cM
        assert_eq!(gm.base_pos(chrom, 2.5), 2_500_000);
        assert_eq!(gm.index2base_pos(chrom, 1), 2_000_000);
    }
}
