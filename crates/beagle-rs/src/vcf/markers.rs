//! Port of `vcf/Markers.java` — an immutable list of markers in chromosome order, with
//! cumulative allele/haplotype-bit sums and bit-packing of haplotype alleles.

use crate::beagleutil::ChromIds;
use crate::blbutil::{consts, BitArray};
use std::collections::HashSet;
use std::hash::{Hash, Hasher};

use super::{marker_utils, Marker};

/// Port of `vcf/Markers.java`.
#[derive(Clone, Debug)]
pub struct Markers {
    markers: Vec<Marker>,
    marker_set: HashSet<Marker>,
    sum_alleles: Vec<i32>,
    sum_hap_bits: Vec<i32>,
}

/// `Integer.SIZE - numberOfLeadingZeros(nAlleles - 1)` — bits to store an allele.
fn n_storage_bits(n_alleles: i32) -> i32 {
    (32 - ((n_alleles - 1) as u32).leading_zeros()) as i32
}

fn check_marker_pos_order(markers: &[Marker]) {
    if markers.len() < 2 {
        return;
    }
    let mut chrom_indices: HashSet<i32> = HashSet::new();
    chrom_indices.insert(markers[0].chrom_index());
    chrom_indices.insert(markers[1].chrom_index());
    for j in 2..markers.len() {
        let chr0 = markers[j - 2].chrom_index();
        let chr1 = markers[j - 1].chrom_index();
        let chr2 = markers[j].chrom_index();
        if chr0 == chr1 && chr1 == chr2 {
            let pos0 = markers[j - 2].pos();
            let pos1 = markers[j - 1].pos();
            let pos2 = markers[j].pos();
            if (pos1 < pos0 && pos1 < pos2) || (pos1 > pos0 && pos1 > pos2) {
                panic!(
                    "markers not in chromosomal order: {}{}{}{}{}{}",
                    consts::NL,
                    markers[j - 2],
                    consts::NL,
                    markers[j - 1],
                    consts::NL,
                    markers[j]
                );
            }
        } else if chr1 != chr2 {
            if chrom_indices.contains(&chr2) {
                panic!(
                    "markers on chromosome are not contiguous: {}",
                    ChromIds::instance().id(chr2)
                );
            }
            chrom_indices.insert(chr2);
        }
    }
}

impl Markers {
    /// `Markers.create(Marker[])` / the private constructor.
    pub fn create(markers: Vec<Marker>) -> Markers {
        check_marker_pos_order(&markers);
        let mut marker_set = HashSet::with_capacity(markers.len());
        for m in &markers {
            if !marker_set.insert(m.clone()) {
                panic!("Duplicate marker: {m}");
            }
        }
        let n = markers.len();
        let mut sum_alleles = vec![0i32; n + 1];
        let mut sum_hap_bits = vec![0i32; n + 1];
        for j in 1..=n {
            sum_alleles[j] = sum_alleles[j - 1] + markers[j - 1].n_alleles();
            sum_hap_bits[j] = sum_hap_bits[j - 1] + n_storage_bits(markers[j - 1].n_alleles());
        }
        Markers {
            markers,
            marker_set,
            sum_alleles,
            sum_hap_bits,
        }
    }

    /// `cumSumGenotypes()`.
    pub fn cum_sum_genotypes(&self) -> Vec<i32> {
        let n = self.markers.len();
        let mut ia = vec![0i32; n + 1];
        for j in 1..=n {
            ia[j] = ia[j - 1] + marker_utils::n_genotypes(self.markers[j - 1].n_alleles());
        }
        ia
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.markers.len() as i32
    }

    /// `marker(int)`.
    pub fn marker(&self, marker: i32) -> &Marker {
        &self.markers[marker as usize]
    }

    /// `markers()`.
    pub fn markers(&self) -> Vec<Marker> {
        self.markers.clone()
    }

    /// `contains(Marker)`.
    pub fn contains(&self, marker: &Marker) -> bool {
        self.marker_set.contains(marker)
    }

    /// `restrict(int start, int end)`.
    pub fn restrict(&self, start: i32, end: i32) -> Markers {
        if end > self.markers.len() as i32 {
            panic!("end > this.nMarkers(): {end}");
        }
        Markers::create(self.markers[start as usize..end as usize].to_vec())
    }

    /// `restrict(int[] indices)` — distinct indices in increasing order.
    pub fn restrict_indices(&self, indices: &[i32]) -> Markers {
        let mut ma = Vec::with_capacity(indices.len());
        ma.push(self.markers[indices[0] as usize].clone());
        for j in 1..indices.len() {
            if indices[j] <= indices[j - 1] {
                panic!("{}", indices[j]);
            }
            ma.push(self.markers[indices[j] as usize].clone());
        }
        Markers::create(ma)
    }

    /// `sumAlleles(int marker)`.
    pub fn sum_alleles(&self, marker: i32) -> i32 {
        self.sum_alleles[marker as usize]
    }

    /// `sumAlleles()`.
    pub fn sum_alleles_total(&self) -> i32 {
        self.sum_alleles[self.markers.len()]
    }

    /// `sumHapBits(int marker)`.
    pub fn sum_hap_bits(&self, marker: i32) -> i32 {
        self.sum_hap_bits[marker as usize]
    }

    /// `sumHapBits()`.
    pub fn sum_hap_bits_total(&self) -> i32 {
        self.sum_hap_bits[self.markers.len()]
    }

    /// `bitsToAlleles(BitArray)`.
    #[allow(clippy::needless_range_loop)] // index spans alleles[m] + sum_hap_bits[m..m+2]
    pub fn bits_to_alleles(&self, hap_bits: &BitArray) -> Vec<i32> {
        let mut alleles = vec![0i32; self.markers.len()];
        for m in 0..self.markers.len() {
            let start = self.sum_hap_bits[m];
            let end = self.sum_hap_bits[m + 1];
            if end == start + 1 {
                alleles[m] = i32::from(hap_bits.get(start));
            } else {
                let mut allele = 0i32;
                let mut mask = 1i32;
                for j in start..end {
                    if hap_bits.get(j) {
                        allele |= mask;
                    }
                    mask <<= 1;
                }
                alleles[m] = allele;
            }
        }
        alleles
    }

    /// `allelesToBits(int[] alleles, BitArray)`.
    #[allow(clippy::needless_range_loop)] // index spans alleles[k] + sum_hap_bits[k..k+2]
    pub fn alleles_to_bits(&self, alleles: &[i32], bit_list: &mut BitArray) {
        assert!(alleles.len() == self.markers.len(), "{}", alleles.len());
        assert!(
            bit_list.size() == self.sum_hap_bits_total(),
            "{}",
            bit_list.size()
        );
        for k in 0..alleles.len() {
            let allele = alleles[k];
            if allele < 0 || allele >= self.markers[k].n_alleles() {
                panic!(
                    "allele \"{}\" out of bounds for marker: {}",
                    allele, self.markers[k]
                );
            }
            let mut mask = 1i32;
            for j in self.sum_hap_bits[k]..self.sum_hap_bits[k + 1] {
                if allele & mask == mask {
                    bit_list.set(j);
                } else {
                    bit_list.clear(j);
                }
                mask <<= 1;
            }
        }
    }

    /// `setAllele(int marker, int allele, BitArray)`.
    pub fn set_allele(&self, marker: i32, allele: i32, bit_list: &mut BitArray) {
        assert!(
            bit_list.size() == self.sum_hap_bits_total(),
            "{}",
            bit_list.size()
        );
        assert!(
            marker >= 0 && (marker as usize) < self.markers.len(),
            "{}",
            marker
        );
        assert!(
            allele >= 0 && allele < self.markers[marker as usize].n_alleles(),
            "{}",
            allele
        );
        let mut mask = 1i32;
        for j in self.sum_hap_bits[marker as usize]..self.sum_hap_bits[marker as usize + 1] {
            if allele & mask == mask {
                bit_list.set(j);
            } else {
                bit_list.clear(j);
            }
            mask <<= 1;
        }
    }

    /// `allele(BitArray hapBits, int marker)`.
    pub fn allele(&self, hap_bits: &BitArray, marker: i32) -> i32 {
        let start = self.sum_hap_bits[marker as usize];
        let end = self.sum_hap_bits[marker as usize + 1];
        if end == start + 1 {
            return i32::from(hap_bits.get(start));
        }
        let mut allele = 0i32;
        let mut mask = 1i32;
        for j in start..end {
            if hap_bits.get(j) {
                allele |= mask;
            }
            mask <<= 1;
        }
        allele
    }
}

impl PartialEq for Markers {
    fn eq(&self, other: &Self) -> bool {
        self.markers == other.markers
    }
}
impl Eq for Markers {}
impl Hash for Markers {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.markers.hash(state);
    }
}

impl std::fmt::Display for Markers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("[")?;
        for (i, m) in self.markers.iter().enumerate() {
            if i > 0 {
                f.write_str(", ")?;
            }
            write!(f, "{m}")?;
        }
        f.write_str("]")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    fn mk(rec: &str) -> Marker {
        Marker::instance(rec, &MarkerParser::new(true, true, true, true))
    }

    fn markers() -> Markers {
        Markers::create(vec![
            mk("chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1"), // 2 alleles -> 1 bit
            mk("chr1\t200\t.\tA\tC,G,T\t.\tPASS\t.\tGT\t0|1"), // 4 alleles -> 2 bits
            mk("chr1\t300\t.\tAC\tA\t.\tPASS\t.\tGT\t0|1"), // 2 alleles -> 1 bit
        ])
    }

    #[test]
    fn sums_and_access() {
        let m = markers();
        assert_eq!(m.size(), 3);
        assert_eq!(m.marker(1).pos(), 200);
        assert!(m.contains(&mk("chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1")));
        assert!(!m.contains(&mk("chr1\t999\t.\tA\tC\t.\tPASS\t.\tGT\t0|1")));
        // sumAlleles: [0, 2, 6, 8]
        assert_eq!(
            (0..=3).map(|k| m.sum_alleles(k)).collect::<Vec<_>>(),
            vec![0, 2, 6, 8]
        );
        // sumHapBits: [0, 1, 3, 4]
        assert_eq!(
            (0..=3).map(|k| m.sum_hap_bits(k)).collect::<Vec<_>>(),
            vec![0, 1, 3, 4]
        );
        assert_eq!(m.cum_sum_genotypes(), vec![0, 3, 13, 16]); // nGeno: 3, 10, 3
    }

    #[test]
    fn bit_packing_roundtrip() {
        let m = markers();
        let mut bits = BitArray::new(m.sum_hap_bits_total());
        let alleles = [1, 3, 0]; // valid for nAlleles 2,4,2
        m.alleles_to_bits(&alleles, &mut bits);
        assert_eq!(m.bits_to_alleles(&bits), alleles.to_vec());
        assert_eq!(m.allele(&bits, 1), 3);
        m.set_allele(1, 2, &mut bits);
        assert_eq!(m.allele(&bits, 1), 2);
    }

    #[test]
    fn restrict_subsets() {
        let m = markers();
        let r = m.restrict(1, 3);
        assert_eq!(r.size(), 2);
        assert_eq!(r.marker(0).pos(), 200);
        let ri = m.restrict_indices(&[0, 2]);
        assert_eq!(ri.size(), 2);
        assert_eq!(ri.marker(1).pos(), 300);
    }
}
