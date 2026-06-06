//! Port of `vcf/Window.java` — a sliding window of target (and optional reference) VCF
//! records, plus the per-allele carrier lists used by imputation.

use std::rc::Rc;

use crate::ints::{IntList, WrappedIntArray};

use super::{BasicGT, GeneticMap, MarkerIndices, RefGT, GT};

/// The classification of an allele's carrier list within a window.
///
/// Java uses two distinct empty `IntArray` sentinels — `Window.ZERO_FREQ_ARRAY` and
/// `Window.HIGH_FREQ_ARRAY` — disambiguated by object identity. Rust models the three
/// states (a real carrier list, no carriers, or too many carriers) as this enum.
#[derive(Clone, Debug)]
pub enum CarrierList {
    /// A list of carrier sample indices (size in `1..=maxCarriers`).
    List(WrappedIntArray),
    /// `Window.ZERO_FREQ_ARRAY` — the allele has no carriers.
    ZeroFreq,
    /// `Window.HIGH_FREQ_ARRAY` — carriers exceed `maxCarriers`.
    HighFreq,
}

/// Port of `vcf/Window.java`.
pub struct Window {
    gen_map: Rc<dyn GeneticMap>,
    window_index: i32,
    last_window: bool,
    indices: MarkerIndices,
    targ_gt: BasicGT,
    ref_gt: Option<RefGT>,
    restrict_ref_gt: Option<RefGT>,
}

impl Window {
    /// `new Window(GeneticMap, int windowIndex, boolean lastWindow, MarkerIndices, RefGT,
    /// BasicGT)`.
    pub fn new(
        gen_map: Rc<dyn GeneticMap>,
        window_index: i32,
        last_window: bool,
        marker_indices: MarkerIndices,
        ref_gt: Option<RefGT>,
        targ_gt: BasicGT,
    ) -> Self {
        assert!(
            targ_gt.n_markers() == marker_indices.n_targ_markers(),
            "inconsistent data"
        );
        if let Some(r) = &ref_gt {
            assert!(
                r.n_markers() == marker_indices.n_markers(),
                "inconsistent data"
            );
        }
        let restrict_ref_gt = ref_gt
            .as_ref()
            .map(|r| r.restrict_to_ref(targ_gt.markers(), &marker_indices.targ_marker_to_marker()));
        Window {
            gen_map,
            window_index,
            last_window,
            indices: marker_indices,
            targ_gt,
            ref_gt,
            restrict_ref_gt,
        }
    }

    /// `genMap()`.
    pub fn gen_map(&self) -> &dyn GeneticMap {
        self.gen_map.as_ref()
    }

    /// A clone of the shared genetic-map handle.
    pub fn gen_map_rc(&self) -> Rc<dyn GeneticMap> {
        self.gen_map.clone()
    }

    /// `chromIndex()`.
    pub fn chrom_index(&self) -> i32 {
        self.targ_gt.marker(0).chrom_index()
    }

    /// `lastWindow()`.
    pub fn last_window(&self) -> bool {
        self.last_window
    }

    /// `windowIndex()` (the first window has index 1).
    pub fn window_index(&self) -> i32 {
        self.window_index
    }

    /// `targGT()`.
    pub fn targ_gt(&self) -> &BasicGT {
        &self.targ_gt
    }

    /// `refGT()`.
    pub fn ref_gt(&self) -> Option<&RefGT> {
        self.ref_gt.as_ref()
    }

    /// `restrictRefGT()`.
    pub fn restrict_ref_gt(&self) -> Option<&RefGT> {
        self.restrict_ref_gt.as_ref()
    }

    /// `indices()`.
    pub fn indices(&self) -> &MarkerIndices {
        &self.indices
    }

    /// `carriers(int maxCarriers)`.
    pub fn carriers(&self, max_carriers: i32) -> Vec<Vec<CarrierList>> {
        (0..self.targ_gt.n_markers())
            .map(|m| self.carriers_at(m, max_carriers))
            .collect()
    }

    fn carriers_at(&self, m: i32, max_carriers: i32) -> Vec<CarrierList> {
        let n_alleles = self.targ_gt.marker(m).n_alleles();
        let mut carriers: Vec<IntList> =
            (0..n_alleles).map(|_| IntList::with_capacity(16)).collect();
        let n_targ_samples = self.targ_gt.n_samples();
        let add = |carriers: &mut [IntList], a1: i32, a2: i32, s: i32| {
            if a1 >= 0 && carriers[a1 as usize].size() <= max_carriers {
                carriers[a1 as usize].add(s);
            }
            if a2 >= 0 && a2 != a1 && carriers[a2 as usize].size() <= max_carriers {
                carriers[a2 as usize].add(s);
            }
        };
        for s in 0..n_targ_samples {
            let hap1 = s << 1;
            let a1 = self.targ_gt.allele(m, hap1);
            let a2 = self.targ_gt.allele(m, hap1 | 0b1);
            add(&mut carriers, a1, a2, s);
        }
        if let Some(ref_gt) = &self.restrict_ref_gt {
            for s in 0..ref_gt.n_samples() {
                let hap1 = s << 1;
                let a1 = ref_gt.allele(m, hap1);
                let a2 = ref_gt.allele(m, hap1 | 0b1);
                add(&mut carriers, a1, a2, n_targ_samples + s);
            }
        }
        carriers
            .into_iter()
            .map(|list| {
                if list.is_empty() {
                    CarrierList::ZeroFreq
                } else if list.size() <= max_carriers {
                    CarrierList::List(WrappedIntArray::from_list(&list))
                } else {
                    CarrierList::HighFreq
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::IntArray;
    use crate::vcf::{
        BasicGTRec, GTRec, MarkerParser, PositionMap, VcfHeader, VcfRecGTParser, HEADER_PREFIX,
    };

    fn header(n_dip: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n_dip {
            hdr.push_str(&format!("\tS{s}"));
        }
        let lines = vec!["##fileformat=VCFv4.2".to_string(), hdr];
        VcfHeader::new_accept_all("src", &lines, &vec![true; n_dip])
    }

    fn rec(h: &VcfHeader, line: &str) -> Rc<dyn GTRec> {
        let p = VcfRecGTParser::new(h, line, &MarkerParser::new(true, true, true, true));
        Rc::new(BasicGTRec::from_parser(&p))
    }

    fn window(lines: &[&str], n_dip: usize) -> Window {
        let h = header(n_dip);
        let targ = BasicGT::new(lines.iter().map(|l| rec(&h, l)).collect());
        let n_targ = targ.n_markers();
        let indices = MarkerIndices::from_counts(0, n_targ, n_targ);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        Window::new(gen_map, 1, true, indices, None, targ)
    }

    #[test]
    fn accessors_and_no_ref() {
        let w = window(
            &[
                "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0",
                "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t0|1",
            ],
            2,
        );
        assert_eq!(w.window_index(), 1);
        assert!(w.last_window());
        assert_eq!(w.targ_gt().n_markers(), 2);
        assert!(w.ref_gt().is_none());
        assert!(w.restrict_ref_gt().is_none());
        assert_eq!(w.chrom_index(), w.targ_gt().marker(0).chrom_index());
    }

    #[test]
    fn carriers_classifies_alleles() {
        // marker0: S0=0|1, S1=1|0 -> allele1 carriers = {0,1}; allele0 carriers = {0,1}
        // marker1: S0=0|0, S1=0|1 -> allele1 carriers = {1}; allele0 carriers = {0,1}
        let w = window(
            &[
                "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0",
                "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t0|1",
            ],
            2,
        );
        let c = w.carriers(10);
        assert_eq!(c.len(), 2);
        // marker0 allele1 -> samples {0,1}
        match &c[0][1] {
            CarrierList::List(a) => assert_eq!(
                (0..a.size()).map(|i| a.get(i)).collect::<Vec<_>>(),
                vec![0, 1]
            ),
            _ => panic!("expected List"),
        }
        // marker1 allele1 -> sample {1}
        match &c[1][1] {
            CarrierList::List(a) => assert_eq!(a.get(0), 1),
            _ => panic!("expected List"),
        }
    }

    #[test]
    fn carriers_high_freq_when_over_max() {
        // 3 samples all carrying allele1 on both haps -> allele1 has 3 carriers
        let w = window(&["chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t1|1\t1|1\t1|1"], 3);
        // maxCarriers=2 -> allele1 (3 carriers) is HighFreq; allele0 (0 carriers) ZeroFreq
        let c = w.carriers(2);
        assert!(matches!(c[0][1], CarrierList::HighFreq));
        assert!(matches!(c[0][0], CarrierList::ZeroFreq));
    }
}
