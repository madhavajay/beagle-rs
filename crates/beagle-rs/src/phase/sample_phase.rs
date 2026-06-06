//! Port of `phase/SamplePhase.java` — an estimated haplotype pair for one sample, with the
//! genotype-cluster structure (homozygous runs split at hets/missing/unphased markers).
//!
//! The static `toBitLists(EstPhase)` factory depends on `EstPhase` and is added with it.

use crate::blbutil::{BitArray, DoubleArray};
use crate::ints::{IntArray, IntList};
use crate::vcf::Markers;

/// Port of `phase/SamplePhase.ClustType`. Discriminants match the Java enum ordinals (the
/// values stored in the `clustType` byte array).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ClustType {
    /// `MISSING_GT`.
    MissingGt = 0,
    /// `MASKED_HET`.
    MaskedHet = 1,
    /// `HOMOZYGOUS_GT`.
    HomozygousGt = 2,
    /// `PHASED_HET`.
    PhasedHet = 3,
    /// `UNPHASED_HET`.
    UnphasedHet = 4,
}

const N_CLUST_TYPES: usize = 5;

impl ClustType {
    /// `ordinal()`.
    pub fn ordinal(self) -> i32 {
        self as i32
    }

    fn from_ordinal(o: i32) -> ClustType {
        match o {
            0 => ClustType::MissingGt,
            1 => ClustType::MaskedHet,
            2 => ClustType::HomozygousGt,
            3 => ClustType::PhasedHet,
            4 => ClustType::UnphasedHet,
            _ => panic!("invalid ClustType ordinal: {o}"),
        }
    }
}

/// Port of `phase/SamplePhase.java`.
pub struct SamplePhase {
    sample: i32,
    markers: Markers,
    hap1: BitArray,
    hap2: BitArray,
    clust_size: Vec<u8>,
    clust_type: Vec<u8>,
    clust_type_cnt: [i32; N_CLUST_TYPES],
}

fn check_increasing(ia: &dyn IntArray, n_markers: i32) {
    let mut last = -1;
    for j in 0..ia.size() {
        assert!(ia.get(j) > last, "non-increasing index list");
        last = ia.get(j);
    }
    assert!(last < n_markers, "index list out of range");
}

fn clust_type_of(is_missing: bool, is_unphased: bool, a1: i32, a2: i32) -> ClustType {
    if is_missing {
        ClustType::MissingGt
    } else if a1 == a2 {
        ClustType::HomozygousGt
    } else if is_unphased {
        ClustType::UnphasedHet
    } else {
        ClustType::PhasedHet
    }
}

#[allow(clippy::too_many_arguments)]
fn set_clusters(
    hap1: &[i32],
    hap2: &[i32],
    missing_gt: &dyn IntArray,
    unph_hets: &dyn IntArray,
    gen_pos: &DoubleArray,
    max_cm: f32,
    clust_type: &mut IntList,
    clust_type_cnt: &mut [i32; N_CLUST_TYPES],
    clust_size_list: &mut IntList,
) {
    let n_markers = gen_pos.size();
    let mut max_clust_end = gen_pos.get(0) + max_cm as f64;
    let mut prev_is_missing_or_het = false;
    let mut last_end = 0;
    let mut miss_index = 0;
    let mut unph_index = 0;
    let mut next_miss = if miss_index < missing_gt.size() {
        let v = missing_gt.get(miss_index);
        miss_index += 1;
        v
    } else {
        -1
    };
    let mut next_unph = if unph_index < unph_hets.size() {
        let v = unph_hets.get(unph_index);
        unph_index += 1;
        v
    } else {
        -1
    };
    let mut prev_type = ClustType::HomozygousGt;
    for m in 0..n_markers {
        let size = m - last_end;
        let ty = clust_type_of(
            m == next_miss,
            m == next_unph,
            hap1[m as usize],
            hap2[m as usize],
        );
        if ty == ClustType::MissingGt {
            next_miss = if miss_index < missing_gt.size() {
                let v = missing_gt.get(miss_index);
                miss_index += 1;
                v
            } else {
                -1
            };
        } else if ty == ClustType::UnphasedHet {
            next_unph = if unph_index < unph_hets.size() {
                let v = unph_hets.get(unph_index);
                unph_index += 1;
                v
            } else {
                -1
            };
        }
        let is_missing_or_het = ty == ClustType::MissingGt
            || ty == ClustType::UnphasedHet
            || ty == ClustType::PhasedHet;
        if is_missing_or_het
            || prev_is_missing_or_het
            || gen_pos.get(m) > max_clust_end
            || size == 255
        {
            if m > 0 {
                clust_type.add(prev_type.ordinal());
                clust_type_cnt[prev_type.ordinal() as usize] += 1;
                clust_size_list.add(size);
                max_clust_end = gen_pos.get(m) + max_cm as f64;
                last_end = m;
            }
            prev_type = ty;
        }
        prev_is_missing_or_het = is_missing_or_het;
    }
    clust_type.add(prev_type.ordinal());
    clust_type_cnt[prev_type.ordinal() as usize] += 1;
    clust_size_list.add(n_markers - last_end);
}

fn to_byte_array(int_list: &IntList) -> Vec<u8> {
    (0..int_list.size())
        .map(|j| int_list.get(j) as u8)
        .collect()
}

impl SamplePhase {
    /// `new SamplePhase(int sample, Markers, DoubleArray genPos, int[] hap1, int[] hap2,
    /// IntArray unphasedHets, IntArray missingGTs)`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        sample: i32,
        markers: Markers,
        gen_pos: &DoubleArray,
        hap1: &[i32],
        hap2: &[i32],
        unphased_hets: &dyn IntArray,
        missing_gts: &dyn IntArray,
    ) -> Self {
        assert!(sample >= 0, "{sample}");
        let n_markers = markers.size();
        assert!(n_markers == gen_pos.size(), "{}", gen_pos.size());
        assert!(hap1.len() as i32 == n_markers, "{}", hap1.len());
        assert!(hap2.len() as i32 == n_markers, "{}", hap2.len());
        check_increasing(unphased_hets, n_markers);
        check_increasing(missing_gts, n_markers);
        let mut bit_hap1 = BitArray::new(markers.sum_hap_bits_total());
        let mut bit_hap2 = BitArray::new(markers.sum_hap_bits_total());
        markers.alleles_to_bits(hap1, &mut bit_hap1);
        markers.alleles_to_bits(hap2, &mut bit_hap2);
        let max_cluster_cm = 0.005f32;
        let mut clust_type_list = IntList::new();
        let mut clust_size_list = IntList::new();
        let mut clust_type_cnt = [0i32; N_CLUST_TYPES];
        set_clusters(
            hap1,
            hap2,
            missing_gts,
            unphased_hets,
            gen_pos,
            max_cluster_cm,
            &mut clust_type_list,
            &mut clust_type_cnt,
            &mut clust_size_list,
        );
        let clust_type = to_byte_array(&clust_type_list);
        let clust_size = to_byte_array(&clust_size_list);
        assert!(clust_size.len() == clust_type.len());
        SamplePhase {
            sample,
            markers,
            hap1: bit_hap1,
            hap2: bit_hap2,
            clust_size,
            clust_type,
            clust_type_cnt,
        }
    }

    /// `sample()`.
    pub fn sample(&self) -> i32 {
        self.sample
    }

    /// `markers()`.
    pub fn markers(&self) -> &Markers {
        &self.markers
    }

    /// `nClusters()`.
    pub fn n_clusters(&self) -> i32 {
        self.clust_size.len() as i32
    }

    /// `clustSize(int cluster)`.
    pub fn clust_size(&self, cluster: i32) -> i32 {
        self.clust_size[cluster as usize] as i32
    }

    /// `clustType(int cluster)`.
    pub fn clust_type(&self, cluster: i32) -> ClustType {
        ClustType::from_ordinal(self.clust_type[cluster as usize] as i32)
    }

    /// `clustEnds()` — exclusive end marker index of each cluster.
    pub fn clust_ends(&self) -> Vec<i32> {
        let mut cum_sum = 0;
        self.clust_size
            .iter()
            .map(|&s| {
                cum_sum += s as i32;
                cum_sum
            })
            .collect()
    }

    /// `nUnphased()`.
    pub fn n_unphased(&self) -> i32 {
        self.clust_type_cnt[ClustType::UnphasedHet.ordinal() as usize]
    }
    /// `nPhased()`.
    pub fn n_phased(&self) -> i32 {
        self.clust_type_cnt[ClustType::PhasedHet.ordinal() as usize]
    }
    /// `nMasked()`.
    pub fn n_masked(&self) -> i32 {
        self.clust_type_cnt[ClustType::MaskedHet.ordinal() as usize]
    }
    /// `nMissing()`.
    pub fn n_missing(&self) -> i32 {
        self.clust_type_cnt[ClustType::MissingGt.ordinal() as usize]
    }
    /// `nHomClusters()`.
    pub fn n_hom_clusters(&self) -> i32 {
        self.clust_type_cnt[ClustType::HomozygousGt.ordinal() as usize]
    }

    /// `maskHetCluster(int cluster)`.
    pub fn mask_het_cluster(&mut self, cluster: i32) {
        assert!(
            self.clust_type[cluster as usize] as i32 == ClustType::UnphasedHet.ordinal(),
            "{}",
            self.clust_type[cluster as usize]
        );
        self.clust_type[cluster as usize] = ClustType::MaskedHet.ordinal() as u8;
        self.clust_type_cnt[ClustType::UnphasedHet.ordinal() as usize] -= 1;
        self.clust_type_cnt[ClustType::MaskedHet.ordinal() as usize] += 1;
    }

    /// `markUnphasedHetClusterAsPhased(int cluster)`.
    pub fn mark_unphased_het_cluster_as_phased(&mut self, cluster: i32) {
        assert!(
            self.clust_type[cluster as usize] as i32 == ClustType::UnphasedHet.ordinal(),
            "{}",
            self.clust_type[cluster as usize]
        );
        self.clust_type[cluster as usize] = ClustType::PhasedHet.ordinal() as u8;
        self.clust_type_cnt[ClustType::UnphasedHet.ordinal() as usize] -= 1;
        self.clust_type_cnt[ClustType::PhasedHet.ordinal() as usize] += 1;
    }

    /// `markMaskedHetClusterAsPhased(int cluster)`.
    pub fn mark_masked_het_cluster_as_phased(&mut self, cluster: i32) {
        assert!(
            self.clust_type[cluster as usize] as i32 == ClustType::MaskedHet.ordinal(),
            "{}",
            self.clust_type[cluster as usize]
        );
        self.clust_type[cluster as usize] = ClustType::PhasedHet.ordinal() as u8;
        self.clust_type_cnt[ClustType::MaskedHet.ordinal() as usize] -= 1;
        self.clust_type_cnt[ClustType::PhasedHet.ordinal() as usize] += 1;
    }

    /// `maskTrailingUnphasedHets()`.
    pub fn mask_trailing_unphased_hets(&mut self) {
        let max_unph_het_clusters = 3;
        let max_masked_base_pairs = 3000;
        let mut unph_het_markers = IntList::new();
        let mut unph_het_clusters = IntList::new();
        let mut start_marker = 0;
        for c in 0..self.clust_type.len() as i32 {
            let ct = self.clust_type(c);
            if ct == ClustType::PhasedHet {
                if 2 <= unph_het_clusters.size()
                    && unph_het_clusters.size() <= max_unph_het_clusters
                {
                    self.mask_trailing(
                        &unph_het_clusters,
                        &unph_het_markers,
                        max_masked_base_pairs,
                    );
                }
                unph_het_markers.clear();
                unph_het_clusters.clear();
            } else if ct == ClustType::UnphasedHet {
                unph_het_markers.add(start_marker);
                unph_het_clusters.add(c);
            }
            start_marker += self.clust_size[c as usize] as i32;
        }
        if 2 <= unph_het_clusters.size() && unph_het_clusters.size() <= max_unph_het_clusters {
            self.mask_trailing(&unph_het_clusters, &unph_het_markers, max_masked_base_pairs);
        }
        debug_assert_eq!(start_marker, self.markers.size());
    }

    fn mask_trailing(
        &mut self,
        unph_het_clusters: &IntList,
        unph_het_markers: &IntList,
        max_masked_base_pairs: i32,
    ) {
        let last_masked_index = unph_het_clusters.size() - 2;
        if last_masked_index == 0 {
            self.mask_het_cluster(unph_het_clusters.get(last_masked_index));
        } else if last_masked_index > 0 {
            let start_pos = self.markers.marker(unph_het_markers.get(0)).pos();
            let end_pos = self
                .markers
                .marker(unph_het_markers.get(last_masked_index))
                .pos();
            if (end_pos - start_pos) <= max_masked_base_pairs {
                for j in 0..=last_masked_index {
                    self.mask_het_cluster(unph_het_clusters.get(j));
                }
            }
        }
    }

    /// `getHaps(BitArray hap1, BitArray hap2)`.
    pub fn get_haps(&self, hap1: &mut BitArray, hap2: &mut BitArray) {
        let n_bits = self.markers.sum_hap_bits_total();
        assert!(
            hap1.size() == n_bits && hap2.size() == n_bits,
            "inconsistent data"
        );
        hap1.copy_from(&self.hap1, 0, self.hap1.size());
        hap2.copy_from(&self.hap2, 0, self.hap2.size());
    }

    /// `allele1(int marker)`.
    pub fn allele1(&self, marker: i32) -> i32 {
        self.markers.allele(&self.hap1, marker)
    }
    /// `allele2(int marker)`.
    pub fn allele2(&self, marker: i32) -> i32 {
        self.markers.allele(&self.hap2, marker)
    }
    /// `setAllele1(int marker, int allele)`.
    pub fn set_allele1(&mut self, marker: i32, allele: i32) {
        self.markers.set_allele(marker, allele, &mut self.hap1);
    }
    /// `setAllele2(int marker, int allele)`.
    pub fn set_allele2(&mut self, marker: i32, allele: i32) {
        self.markers.set_allele(marker, allele, &mut self.hap2);
    }

    /// `swapHaps(int start, int end)`.
    pub fn swap_haps(&mut self, start: i32, end: i32) {
        let start_bit = self.markers.sum_hap_bits(start);
        let end_bit = self.markers.sum_hap_bits(end);
        BitArray::swap_bits(&mut self.hap1, &mut self.hap2, start_bit, end_bit);
    }

    /// `hap1()` — a copy of the first haplotype.
    pub fn hap1(&self) -> BitArray {
        self.hap1.clone()
    }
    /// `hap2()` — a copy of the second haplotype.
    pub fn hap2(&self) -> BitArray {
        self.hap2.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::WrappedIntArray;
    use crate::vcf::{Marker, MarkerParser};

    fn markers(positions: &[i32]) -> Markers {
        let mp = MarkerParser::new(true, true, true, true);
        let ms: Vec<Marker> = positions
            .iter()
            .map(|p| Marker::instance(&format!("chr1\t{p}\t.\tA\tC\t.\tPASS\t.\tGT\t0|0"), &mp))
            .collect();
        Markers::create(ms)
    }

    fn sample_phase() -> SamplePhase {
        // 3 biallelic markers, close together (one homozygous run except a het in the middle)
        let m = markers(&[100, 200, 300]);
        let gen_pos = DoubleArray::from_doubles(&[0.0, 1e-4, 2e-4]);
        let hap1 = [0, 0, 1];
        let hap2 = [0, 1, 1];
        let empty = WrappedIntArray::from_slice(&[]);
        SamplePhase::new(0, m, &gen_pos, &hap1, &hap2, &empty, &empty)
    }

    #[test]
    fn alleles_and_clusters() {
        let sp = sample_phase();
        assert_eq!(sp.sample(), 0);
        assert_eq!(sp.allele1(0), 0);
        assert_eq!(sp.allele1(2), 1);
        assert_eq!(sp.allele2(1), 1);
        // marker1 is a phased het -> splits the homozygous run into 3 clusters
        assert_eq!(sp.n_clusters(), 3);
        assert_eq!(sp.clust_size(0) + sp.clust_size(1) + sp.clust_size(2), 3);
        assert_eq!(sp.clust_ends(), vec![1, 2, 3]);
        assert_eq!(sp.clust_type(1), ClustType::PhasedHet);
        assert_eq!(sp.n_phased(), 1);
    }

    #[test]
    fn set_allele_and_swap() {
        let mut sp = sample_phase();
        sp.set_allele1(0, 1);
        assert_eq!(sp.allele1(0), 1);
        // swap all markers: hap1<->hap2
        let before1: Vec<i32> = (0..3).map(|m| sp.allele1(m)).collect();
        let before2: Vec<i32> = (0..3).map(|m| sp.allele2(m)).collect();
        sp.swap_haps(0, 3);
        let after1: Vec<i32> = (0..3).map(|m| sp.allele1(m)).collect();
        let after2: Vec<i32> = (0..3).map(|m| sp.allele2(m)).collect();
        assert_eq!(after1, before2);
        assert_eq!(after2, before1);
    }

    #[test]
    fn hap_copies_are_independent() {
        let sp = sample_phase();
        let h1 = sp.hap1();
        assert_eq!(h1.size(), sp.markers().sum_hap_bits_total());
    }
}
