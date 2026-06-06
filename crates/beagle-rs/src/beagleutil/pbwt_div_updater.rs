//! Port of `beagleutil/PbwtDivUpdater.java` — updates PBWT prefix AND divergence
//! arrays (forward and backward). Durbin (2014), Bioinformatics 30(9):1266-1272.

use crate::ints::{IntArray, IntList};

/// Port of `beagleutil/PbwtDivUpdater.java`.
pub struct PbwtDivUpdater {
    n_haps: i32,
    a: Vec<IntList>, // prefix buckets per allele
    d: Vec<IntList>, // divergence buckets per allele
    p: Vec<i32>,
}

const INIT_NUM_ALLELES: usize = 4;

impl PbwtDivUpdater {
    /// `new PbwtDivUpdater(int nHaps)`.
    pub fn new(n_haps: i32) -> Self {
        PbwtDivUpdater {
            n_haps,
            a: (0..INIT_NUM_ALLELES).map(|_| IntList::new()).collect(),
            d: (0..INIT_NUM_ALLELES).map(|_| IntList::new()).collect(),
            p: vec![0; INIT_NUM_ALLELES],
        }
    }

    /// `nHaps()`.
    pub fn n_haps(&self) -> i32 {
        self.n_haps
    }

    /// `fwdUpdate(rec, nAlleles, marker, prefix, div)` — forward PBWT update.
    pub fn fwd_update(
        &mut self,
        rec: &dyn IntArray,
        n_alleles: i32,
        marker: i32,
        prefix: &mut [i32],
        div: &mut [i32],
    ) {
        assert!(rec.size() == self.n_haps, "{}", rec.size());
        assert!(prefix.len() as i32 == self.n_haps, "{}", prefix.len());
        self.init_p_array(n_alleles, marker + 1);
        for i in 0..self.n_haps as usize {
            let allele = rec.get(prefix[i]);
            assert!(allele < n_alleles, "{}", n_alleles);
            let allele = allele as usize;
            for j in 0..n_alleles as usize {
                if div[i] > self.p[j] {
                    self.p[j] = div[i];
                }
            }
            self.a[allele].add(prefix[i]);
            self.d[allele].add(self.p[allele]);
            self.p[allele] = i32::MIN;
        }
        self.update_prefix_and_div(n_alleles, prefix, div);
    }

    /// `bwdUpdate(rec, nAlleles, marker, prefix, div)` — backward PBWT update.
    pub fn bwd_update(
        &mut self,
        rec: &dyn IntArray,
        n_alleles: i32,
        marker: i32,
        prefix: &mut [i32],
        div: &mut [i32],
    ) {
        assert!(rec.size() == self.n_haps, "{}", rec.size());
        assert!(prefix.len() as i32 == self.n_haps, "{}", prefix.len());
        self.init_p_array(n_alleles, marker - 1);
        for i in 0..self.n_haps as usize {
            let allele = rec.get(prefix[i]);
            assert!(allele < n_alleles, "{}", n_alleles);
            let allele = allele as usize;
            for j in 0..n_alleles as usize {
                if div[i] < self.p[j] {
                    self.p[j] = div[i];
                }
            }
            self.a[allele].add(prefix[i]);
            self.d[allele].add(self.p[allele]);
            self.p[allele] = i32::MAX;
        }
        self.update_prefix_and_div(n_alleles, prefix, div);
    }

    fn update_prefix_and_div(&mut self, n_alleles: i32, prefix: &mut [i32], div: &mut [i32]) {
        let mut start = 0usize;
        for al in 0..n_alleles as usize {
            let pa = self.a[al].to_array();
            let da = self.d[al].to_array();
            prefix[start..start + pa.len()].copy_from_slice(&pa);
            div[start..start + da.len()].copy_from_slice(&da);
            start += pa.len();
            self.a[al].clear();
            self.d[al].clear();
        }
        debug_assert_eq!(start as i32, self.n_haps);
    }

    fn init_p_array(&mut self, n_alleles: i32, init_value: i32) {
        assert!(n_alleles >= 1, "{}", n_alleles);
        if n_alleles as usize > self.a.len() {
            self.a.resize_with(n_alleles as usize, IntList::new);
            self.d.resize_with(n_alleles as usize, IntList::new);
            self.p.resize(n_alleles as usize, 0);
        }
        for v in self.p[0..n_alleles as usize].iter_mut() {
            *v = init_value;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::WrappedIntArray;

    #[test]
    fn forward_prefix_and_div_update() {
        let rec = WrappedIntArray::from_slice(&[1, 0, 1, 0]);
        let mut u = PbwtDivUpdater::new(4);
        let mut prefix = [0, 1, 2, 3];
        let mut div = [0, 0, 0, 0];
        u.fwd_update(&rec, 2, 5, &mut prefix, &mut div);
        // traced by hand against the Java algorithm
        assert_eq!(prefix, [1, 3, 0, 2]);
        assert_eq!(div, [6, 0, 6, 0]);
    }
}
