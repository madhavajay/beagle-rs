//! Port of `main/Pedigree.java` — parent-offspring relationships within a sample list
//! (singles, duos, trios), read from an optional linkage-format pedigree file.

use std::path::Path;

use crate::beagleutil::SampleIds;
use crate::blbutil::{InputIt, StringUtil};
use crate::vcf::Samples;

const NO_PARENT: &str = "0";

/// Port of `main/Pedigree.java`.
#[derive(Clone)]
pub struct Pedigree {
    samples: Samples,
    singles: Vec<i32>,
    relateds: Vec<i32>,
    duo_offspring: Vec<i32>,
    trio_offspring: Vec<i32>,
    mothers: Vec<i32>,
    fathers: Vec<i32>,
    offspring: Vec<Vec<i32>>,
}

fn n_parents(index: usize, father: &[i32], mother: &[i32]) -> i32 {
    let mut cnt = 0;
    if father[index] >= 0 {
        cnt += 1;
    }
    if mother[index] >= 0 {
        cnt += 1;
    }
    cnt
}

fn index(id_index_to_index: &[i32], id: &str) -> i32 {
    let id_index = SampleIds::instance().get_index_if_indexed(id);
    if id_index != -1 && (id_index as usize) < id_index_to_index.len() {
        id_index_to_index[id_index as usize]
    } else {
        -1
    }
}

fn get_ped_fields(line: &str) -> Vec<&str> {
    let fields = StringUtil::get_fields_ws_limit(line, 5);
    assert!(fields.len() >= 4, "invalid line in ped file: {line}");
    fields
}

fn set_parent_offspring(
    parent: i32,
    child: i32,
    sample2_parent: &mut [i32],
    children: &mut [Option<Vec<i32>>],
) {
    if parent != -1 {
        sample2_parent[child as usize] = parent;
        children[parent as usize]
            .get_or_insert_with(|| Vec::with_capacity(3))
            .push(child);
    }
}

fn read_ped_line(
    id_index_to_index: &[i32],
    processed: &mut [bool],
    line: &str,
    fathers: &mut [i32],
    mothers: &mut [i32],
    children: &mut [Option<Vec<i32>>],
) {
    if !line.is_empty() {
        let sa = get_ped_fields(line);
        let child_id = sa[1];
        let child = index(id_index_to_index, child_id);
        if child != -1 {
            assert!(
                !processed[child as usize],
                "duplicate sample in pedigree file: {child_id}"
            );
            processed[child as usize] = true;
            let father = if sa[2] == NO_PARENT {
                -1
            } else {
                index(id_index_to_index, sa[2])
            };
            let mother = if sa[3] == NO_PARENT {
                -1
            } else {
                index(id_index_to_index, sa[3])
            };
            set_parent_offspring(father, child, fathers, children);
            set_parent_offspring(mother, child, mothers, children);
        }
    }
}

fn read_ped_file(
    samples: &Samples,
    ped_file: Option<&Path>,
    fathers: &mut [i32],
    mothers: &mut [i32],
    offspring: &mut [Vec<i32>],
) {
    if let Some(ped_file) = ped_file {
        let id_index_to_index = SampleIds::instance().id_index_to_index(&samples.ids());
        let mut processed = vec![false; samples.size() as usize];
        let mut children: Vec<Option<Vec<i32>>> = vec![None; samples.size() as usize];
        for line in InputIt::from_gzip_file(ped_file) {
            let line = line.trim();
            read_ped_line(
                &id_index_to_index,
                &mut processed,
                line,
                fathers,
                mothers,
                &mut children,
            );
        }
        for (j, child) in children.into_iter().enumerate() {
            if let Some(mut ia) = child {
                ia.sort_unstable();
                offspring[j] = ia;
            }
        }
    }
}

impl Pedigree {
    /// `new Pedigree(Samples samples, File pedFile)`.
    pub fn new(samples: Samples, ped_file: Option<&Path>) -> Self {
        let n_samples = samples.size() as usize;
        let mut fathers = vec![-1i32; n_samples];
        let mut mothers = vec![-1i32; n_samples];
        let mut offspring: Vec<Vec<i32>> = vec![Vec::new(); n_samples];

        read_ped_file(
            &samples,
            ped_file,
            &mut fathers,
            &mut mothers,
            &mut offspring,
        );

        let mut singles = Vec::new();
        let mut duo_offspring = Vec::new();
        let mut trio_offspring = Vec::new();
        let mut relateds = Vec::new();
        for (s, offspring_s) in offspring.iter().enumerate() {
            match n_parents(s, &fathers, &mothers) {
                0 => {
                    if offspring_s.is_empty() {
                        singles.push(s as i32);
                    } else {
                        relateds.push(s as i32);
                    }
                }
                1 => {
                    duo_offspring.push(s as i32);
                    relateds.push(s as i32);
                }
                2 => {
                    trio_offspring.push(s as i32);
                    relateds.push(s as i32);
                }
                _ => unreachable!(),
            }
        }

        Pedigree {
            samples,
            singles,
            relateds,
            duo_offspring,
            trio_offspring,
            mothers,
            fathers,
            offspring,
        }
    }

    /// `samples()`.
    pub fn samples(&self) -> &Samples {
        &self.samples
    }
    /// `nSamples()`.
    pub fn n_samples(&self) -> i32 {
        self.samples.size()
    }
    /// `nSingles()`.
    pub fn n_singles(&self) -> i32 {
        self.singles.len() as i32
    }
    /// `nDuos()`.
    pub fn n_duos(&self) -> i32 {
        self.duo_offspring.len() as i32
    }
    /// `nTrios()`.
    pub fn n_trios(&self) -> i32 {
        self.trio_offspring.len() as i32
    }
    /// `singles()`.
    pub fn singles(&self) -> Vec<i32> {
        self.singles.clone()
    }
    /// `relateds()`.
    pub fn relateds(&self) -> Vec<i32> {
        self.relateds.clone()
    }
    /// `single(int index)`.
    pub fn single(&self, index: i32) -> i32 {
        self.singles[index as usize]
    }
    /// `duoParent(int index)`.
    pub fn duo_parent(&self, index: i32) -> i32 {
        let offspring = self.duo_offspring[index as usize] as usize;
        if self.fathers[offspring] >= 0 {
            self.fathers[offspring]
        } else {
            debug_assert!(self.mothers[offspring] >= 0);
            self.mothers[offspring]
        }
    }
    /// `duoOffspring(int index)`.
    pub fn duo_offspring(&self, index: i32) -> i32 {
        self.duo_offspring[index as usize]
    }
    /// `trioFather(int index)`.
    pub fn trio_father(&self, index: i32) -> i32 {
        self.fathers[self.trio_offspring[index as usize] as usize]
    }
    /// `trioMother(int index)`.
    pub fn trio_mother(&self, index: i32) -> i32 {
        self.mothers[self.trio_offspring[index as usize] as usize]
    }
    /// `trioOffspring(int index)`.
    pub fn trio_offspring(&self, index: i32) -> i32 {
        self.trio_offspring[index as usize]
    }
    /// `father(int sample)`.
    pub fn father(&self, sample: i32) -> i32 {
        self.fathers[sample as usize]
    }
    /// `mother(int sample)`.
    pub fn mother(&self, sample: i32) -> i32 {
        self.mothers[sample as usize]
    }
    /// `nOffspring(int sample)`.
    pub fn n_offspring(&self, sample: i32) -> i32 {
        self.offspring[sample as usize].len() as i32
    }
    /// `offspring(int sample, int index)`.
    pub fn offspring(&self, sample: i32, index: i32) -> i32 {
        self.offspring[sample as usize][index as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_ped_file_all_singles() {
        let samples = Samples::new(
            &["NPED0".into(), "NPED1".into(), "NPED2".into()],
            &[true; 3],
        );
        let ped = Pedigree::new(samples, None);
        assert_eq!(ped.n_samples(), 3);
        assert_eq!(ped.n_singles(), 3);
        assert_eq!(ped.n_duos(), 0);
        assert_eq!(ped.n_trios(), 0);
        assert_eq!(ped.singles(), vec![0, 1, 2]);
        assert!(ped.relateds().is_empty());
        assert_eq!(ped.father(0), -1);
        assert_eq!(ped.mother(2), -1);
        assert_eq!(ped.n_offspring(1), 0);
    }

    #[test]
    fn reads_trio_and_duo_from_ped_file() {
        use std::io::Write;
        // register these ids in the global SampleIds via Samples::new
        let samples = Samples::new(
            &[
                "kid".into(),
                "dad".into(),
                "mom".into(),
                "duokid".into(),
                "lone".into(),
            ],
            &[true; 5],
        );
        // famID indivID fatherID motherID
        let ped = "\
F1\tkid\tdad\tmom
F1\tdad\t0\t0
F1\tmom\t0\t0
F2\tduokid\tlone\t0
F2\tlone\t0\t0
";
        let path = std::env::temp_dir().join("beagle_rs_test.ped");
        std::fs::File::create(&path)
            .unwrap()
            .write_all(ped.as_bytes())
            .unwrap();

        let p = Pedigree::new(samples, Some(&path));
        assert_eq!(p.n_trios(), 1);
        assert_eq!(p.n_duos(), 1);
        // trio offspring is "kid" (index 0); father "dad" (1), mother "mom" (2)
        assert_eq!(p.trio_offspring(0), 0);
        assert_eq!(p.trio_father(0), 1);
        assert_eq!(p.trio_mother(0), 2);
        // duo offspring is "duokid" (3); parent "lone" (4)
        assert_eq!(p.duo_offspring(0), 3);
        assert_eq!(p.duo_parent(0), 4);
        // dad/mom/lone have offspring -> relateds, not singles
        assert_eq!(p.n_singles(), 0);
        assert_eq!(p.n_offspring(1), 1); // dad -> kid
        assert_eq!(p.offspring(1, 0), 0);
    }
}
