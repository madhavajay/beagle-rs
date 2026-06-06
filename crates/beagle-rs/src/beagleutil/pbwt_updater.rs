//! Port of `beagleutil/PbwtUpdater.java` — updates a PBWT prefix array (no divergence
//! array). Durbin (2014), Bioinformatics 30(9):1266-1272.

use crate::ints::{IntArray, IntList};

/// Port of `beagleutil/PbwtUpdater.java`.
pub struct PbwtUpdater {
    n_haps: i32,
    a: Vec<IntList>, // one bucket per allele
}

const INIT_NUM_ALLELES: usize = 4;

impl PbwtUpdater {
    /// `new PbwtUpdater(int nHaps)`.
    pub fn new(n_haps: i32) -> Self {
        assert!(n_haps >= 0, "{}", n_haps);
        PbwtUpdater {
            n_haps,
            a: (0..INIT_NUM_ALLELES).map(|_| IntList::new()).collect(),
        }
    }

    /// `nHaps()`.
    pub fn n_haps(&self) -> i32 {
        self.n_haps
    }

    /// `update(IntArray rec, int nAlleles, int[] prefix)` — forward PBWT prefix update.
    pub fn update(&mut self, rec: &dyn IntArray, n_alleles: i32, prefix: &mut [i32]) {
        assert!(rec.size() == self.n_haps, "{}", rec.size());
        assert!(prefix.len() as i32 == self.n_haps, "{}", prefix.len());
        self.initialize_arrays(n_alleles);
        for &h in prefix.iter() {
            let allele = rec.get(h);
            assert!(allele < n_alleles, "{}", allele);
            self.a[allele as usize].add(h);
        }
        self.update_prefix(n_alleles, prefix);
    }

    /// `update(int[] alleles, int nAlleles, int[] prefix)`.
    pub fn update_alleles(&mut self, alleles: &[i32], n_alleles: i32, prefix: &mut [i32]) {
        assert!(alleles.len() as i32 == self.n_haps, "{}", alleles.len());
        assert!(prefix.len() as i32 == self.n_haps, "{}", prefix.len());
        self.initialize_arrays(n_alleles);
        for &h in prefix.iter() {
            let allele = alleles[h as usize];
            assert!(allele < n_alleles, "{}", allele);
            self.a[allele as usize].add(h);
        }
        self.update_prefix(n_alleles, prefix);
    }

    fn update_prefix(&mut self, n_alleles: i32, prefix: &mut [i32]) {
        let mut start = 0usize;
        for al in 0..n_alleles as usize {
            let arr = self.a[al].to_array();
            prefix[start..start + arr.len()].copy_from_slice(&arr);
            start += arr.len();
            self.a[al].clear();
        }
        debug_assert_eq!(start as i32, self.n_haps);
    }

    fn initialize_arrays(&mut self, n_alleles: i32) {
        assert!(n_alleles >= 1, "{}", n_alleles);
        if n_alleles as usize > self.a.len() {
            self.a.resize_with(n_alleles as usize, IntList::new);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::WrappedIntArray;

    #[test]
    fn forward_prefix_update() {
        // haps 0..4 with alleles [1,0,1,0]; stable partition by allele.
        let rec = WrappedIntArray::from_slice(&[1, 0, 1, 0]);
        let mut u = PbwtUpdater::new(4);
        let mut prefix = [0, 1, 2, 3];
        u.update(&rec, 2, &mut prefix);
        // allele-0 haps (1,3) then allele-1 haps (0,2), in prefix order.
        assert_eq!(prefix, [1, 3, 0, 2]);

        // int[] variant agrees
        let mut u2 = PbwtUpdater::new(4);
        let mut prefix2 = [0, 1, 2, 3];
        u2.update_alleles(&[1, 0, 1, 0], 2, &mut prefix2);
        assert_eq!(prefix2, [1, 3, 0, 2]);
    }

    #[test]
    fn second_update_uses_prior_prefix() {
        let rec = WrappedIntArray::from_slice(&[0, 0, 1, 1]);
        let mut u = PbwtUpdater::new(4);
        let mut prefix = [3, 2, 1, 0];
        u.update(&rec, 2, &mut prefix);
        // prefix order [3,2,1,0]: alleles 1,1,0,0 -> allele0 haps (1,0), allele1 (3,2)
        assert_eq!(prefix, [1, 0, 3, 2]);
    }
}
