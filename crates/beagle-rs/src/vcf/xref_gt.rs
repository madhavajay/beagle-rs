//! Port of `vcf/XRefGT.java` — phased, non-missing genotypes stored in haplotype-major
//! (column-major) order: one `BitArray` per haplotype packing every marker's allele bits.
//!
//! The parallelism in Java's `from`/`fromPhasedGT` only partitions work; the output array
//! order is preserved, so the Rust port builds it sequentially (identical result) and
//! ignores the thread count.

use crate::blbutil::{consts, BitArray};
use crate::phase::SamplePhase;
use std::fmt;
use std::rc::Rc;

use super::{Marker, Markers, RestrictedGT, Samples, GT};

/// Port of `vcf/XRefGT.java`.
#[derive(Clone)]
pub struct XRefGT {
    samples: Samples,
    markers: Markers,
    haps: Vec<BitArray>,
}

impl XRefGT {
    /// Private constructor (Java `XRefGT(Markers, Samples, BitArray[])`): callers must
    /// validate arguments and ensure `haps` does not alias.
    fn new_parts(markers: Markers, samples: Samples, haps: Vec<BitArray>) -> Self {
        XRefGT {
            samples,
            markers,
            haps,
        }
    }

    /// `XRefGT.combine(XRefGT first, XRefGT second)` — concatenates the haplotypes of two
    /// instances over the same markers; sample lists must be disjoint.
    pub fn combine(first: &XRefGT, second: &XRefGT) -> XRefGT {
        assert!(first.markers == second.markers, "inconsisent data");
        let samples = Samples::combine(&first.samples, &second.samples);
        let mut haps = Vec::with_capacity(first.haps.len() + second.haps.len());
        haps.extend(first.haps.iter().cloned());
        haps.extend(second.haps.iter().cloned());
        XRefGT::new_parts(first.markers.clone(), samples, haps)
    }

    /// `XRefGT.from(Samples, AtomicReferenceArray<SamplePhase>)` — builds haplotype-major
    /// genotypes from per-sample phasings (`phase[s]` supplies sample `s`'s two haplotypes).
    /// Java parallelizes over haplotypes, but the output order is preserved, so the Rust
    /// port fills the array sequentially.
    pub fn from(samples: &Samples, phase: &[SamplePhase]) -> XRefGT {
        let n_samples = phase.len();
        assert!(
            n_samples != 0 && samples.size() as usize == n_samples,
            "{n_samples}"
        );
        let markers = phase[0].markers().clone();
        let n_haps = n_samples << 1;
        let haps: Vec<BitArray> = (0..n_haps)
            .map(|h| {
                let samp_phase = &phase[h >> 1];
                assert!(*samp_phase.markers() == markers, "inconsistent data");
                if (h & 0b1) == 0 {
                    samp_phase.hap1()
                } else {
                    samp_phase.hap2()
                }
            })
            .collect();
        XRefGT::new_parts(markers, samples.clone(), haps)
    }

    /// `XRefGT.fromPhasedGT(GT gt, int nThreads)`. `n_threads` does not affect output.
    pub fn from_phased_gt(gt: &dyn GT, n_threads: i32) -> XRefGT {
        assert!(n_threads >= 1, "{n_threads}");
        let haps = hap_data(gt);
        XRefGT::new_parts(gt.markers().clone(), gt.samples().clone(), haps)
    }

    /// `hash(int hap, int start, int end)` — hash of a haplotype's alleles over markers
    /// `[start, end)`.
    pub fn hash(&self, hap: i32, start: i32, end: i32) -> i32 {
        let start_bit = self.markers.sum_hap_bits(start);
        let end_bit = self.markers.sum_hap_bits(end);
        self.haps[hap as usize].hash(start_bit, end_bit)
    }

    /// `copyTo(int hap, int start, int end, BitArray bitList)`.
    pub fn copy_to(&self, hap: i32, start: i32, end: i32, bit_list: &mut BitArray) {
        let start_bit = self.markers.sum_hap_bits(start);
        let end_bit = self.markers.sum_hap_bits(end);
        bit_list.copy_from(&self.haps[hap as usize], start_bit, end_bit);
    }

    /// `restrict(int start, int end)` returning the concrete `XRefGT`.
    pub fn restrict_to_xref(&self, start: i32, end: i32) -> XRefGT {
        let restrict_markers = self.markers.restrict(start, end);
        let start_bit = self.markers.sum_hap_bits(start);
        let end_bit = self.markers.sum_hap_bits(end);
        let restrict_haps: Vec<BitArray> = self
            .haps
            .iter()
            .map(|h| h.restrict(start_bit, end_bit))
            .collect();
        XRefGT::new_parts(restrict_markers, self.samples.clone(), restrict_haps)
    }
}

fn hap_data(gt: &dyn GT) -> Vec<BitArray> {
    assert!(gt.is_phased(), "{gt}", gt = DisplayGt(gt));
    let n_hap_bits = gt.markers().sum_hap_bits_total();
    let n_haps = gt.n_haps();
    let n_markers = gt.n_markers();
    let mut haps: Vec<BitArray> = (0..n_haps).map(|_| BitArray::new(n_hap_bits)).collect();
    for m in 0..n_markers {
        set_hap_data(gt, m, &mut haps);
    }
    haps
}

fn set_hap_data(phased_gt: &dyn GT, m: i32, haps: &mut [BitArray]) {
    let markers = phased_gt.markers();
    let start_bit = markers.sum_hap_bits(m);
    let end_bit = markers.sum_hap_bits(m + 1);
    for (j, hap) in haps.iter_mut().enumerate() {
        let allele = phased_gt.allele(m, j as i32);
        // LSB-first packing across [start_bit, end_bit), matching Markers::allele decoding.
        let mut mask = 1i32;
        for i in start_bit..end_bit {
            if (allele & mask) == mask {
                hap.set(i);
            }
            mask <<= 1;
        }
    }
}

/// Helper so `assert!`'s panic message can use `gt`'s `GT` (which is not `Display`).
struct DisplayGt<'a>(&'a dyn GT);
impl fmt::Display for DisplayGt<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "GT(nMarkers={}, nSamples={})",
            self.0.n_markers(),
            self.0.n_samples()
        )
    }
}

impl GT for XRefGT {
    fn is_reversed(&self) -> bool {
        false
    }

    fn n_markers(&self) -> i32 {
        self.markers.size()
    }

    fn marker(&self, marker_index: i32) -> &Marker {
        self.markers.marker(marker_index)
    }

    fn markers(&self) -> &Markers {
        &self.markers
    }

    fn n_haps(&self) -> i32 {
        self.haps.len() as i32
    }

    fn n_samples(&self) -> i32 {
        self.samples.size()
    }

    fn samples(&self) -> &Samples {
        &self.samples
    }

    fn is_phased(&self) -> bool {
        true
    }

    fn allele(&self, marker: i32, hap: i32) -> i32 {
        self.markers.allele(&self.haps[hap as usize], marker)
    }

    fn restrict(self: Rc<Self>, markers: &Markers, indices: &[i32]) -> Rc<dyn GT> {
        Rc::new(RestrictedGT::new(self, markers.clone(), indices))
    }

    fn restrict_range(self: Rc<Self>, start: i32, end: i32) -> Rc<dyn GT> {
        Rc::new(self.restrict_to_xref(start, end))
    }
}

impl fmt::Display for XRefGT {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "class vcf.XRefGT nMarkers={} nSamples={}{}",
            self.n_markers(),
            self.n_samples(),
            consts::NL
        )?;
        let n_markers = self.markers.size();
        let n_samples = self.samples.size();
        for m in 0..n_markers {
            write!(f, "{}", self.markers.marker(m))?;
            write!(
                f,
                "{nl}{miss}{tab}PASS{tab}{miss}{tab}GT",
                nl = consts::NL,
                miss = consts::MISSING_DATA_CHAR,
                tab = consts::TAB
            )?;
            for s in 0..n_samples {
                let hap1 = s << 1;
                write!(
                    f,
                    "{}{}{}{}",
                    consts::TAB,
                    self.allele(m, hap1),
                    consts::PHASED_SEP,
                    self.allele(m, hap1 | 0b1)
                )?;
            }
        }
        write!(f, "{}]", consts::NL)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, MarkerParser, VcfHeader, VcfRecGTParser, HEADER_PREFIX,
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

    fn phased_gt(h: &VcfHeader, lines: &[&str]) -> BasicGT {
        BasicGT::new(lines.iter().map(|l| rec(h, l)).collect())
    }

    #[test]
    fn round_trips_alleles_through_bit_packing() {
        let h = header(2);
        // biallelic + triallelic markers to exercise multi-bit packing
        let gt = phased_gt(
            &h,
            &[
                "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0",
                "chr1\t200\t.\tG\tT,A\t.\tPASS\t.\tGT\t2|1\t0|2",
            ],
        );
        let x = XRefGT::from_phased_gt(&gt, 1);
        assert_eq!(x.n_markers(), 2);
        assert_eq!(x.n_haps(), 4);
        assert!(x.is_phased());
        for m in 0..2 {
            for hp in 0..4 {
                assert_eq!(x.allele(m, hp), gt.allele(m, hp), "m={m} hap={hp}");
            }
        }
    }

    #[test]
    fn restrict_and_combine() {
        let h = header(1);
        let gt = phased_gt(
            &h,
            &[
                "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
                "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|0",
                "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1",
            ],
        );
        let x = XRefGT::from_phased_gt(&gt, 1);
        let sub = x.restrict_to_xref(1, 3);
        assert_eq!(sub.n_markers(), 2);
        assert_eq!(sub.marker(0).pos(), 200);
        assert_eq!(sub.allele(0, 0), 1);
        assert_eq!(sub.allele(1, 1), 1);

        // combine two single-sample XRefGTs over the same markers
        let h2 = header(1);
        let gt2 = phased_gt(
            &h2,
            &[
                "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t1|1",
                "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0",
                "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|1",
            ],
        );
        let x2 = XRefGT::from_phased_gt(&gt2, 1);
        let combined = XRefGT::combine(&x, &x2);
        assert_eq!(combined.n_haps(), 4); // 2 + 2
        assert_eq!(combined.n_samples(), 2);
        // first two haps from x, next two from x2
        assert_eq!(combined.allele(0, 0), 0);
        assert_eq!(combined.allele(0, 2), 1); // x2 sample hap0 at marker0
    }

    #[test]
    fn hash_matches_underlying_bitarray() {
        let h = header(1);
        let gt = phased_gt(
            &h,
            &[
                "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
                "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|0",
            ],
        );
        let x = XRefGT::from_phased_gt(&gt, 1);
        // hash over both markers for hap 0 should equal the BitArray hash of those bits
        let hv = x.hash(0, 0, 2);
        let mut copy = BitArray::new(x.markers.sum_hap_bits_total());
        x.copy_to(0, 0, 2, &mut copy);
        assert_eq!(hv, copy.hash(0, x.markers.sum_hap_bits(2)));
    }
}
