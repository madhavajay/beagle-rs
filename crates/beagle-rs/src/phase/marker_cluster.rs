//! Port of `phase/MarkerCluster.java` — partitions a sample's markers into the contiguous
//! genotype clusters of its current `SamplePhase`, exposing per-cluster boundaries, the
//! unphased-het cluster list, and the per-cluster recombination probability.

use crate::blbutil::FloatArray;
use crate::ints::WrappedIntArray;

use super::{ClustType, PhaseData, SamplePhase};

/// Port of `phase/MarkerCluster.java`.
pub struct MarkerCluster {
    sample_phase: SamplePhase,
    cluster_to_end: Vec<i32>,
    unph_het_clusters: WrappedIntArray,
    p_recomb: FloatArray,
}

fn unph_het_clusters(sample_phase: &SamplePhase) -> WrappedIntArray {
    let n_unph = sample_phase.n_unphased();
    let n_clusters = sample_phase.n_clusters();
    let mut clusters = vec![0i32; n_unph as usize];
    let mut index = 0;
    for c in 0..n_clusters {
        if sample_phase.clust_type(c) == ClustType::UnphasedHet {
            clusters[index] = c;
            index += 1;
        }
    }
    WrappedIntArray::from_slice(&clusters)
}

fn p_clust_recomb(p_recomb: &FloatArray, cluster_to_end: &[i32]) -> FloatArray {
    let n_clusters = cluster_to_end.len();
    let mut p_clust_recomb = vec![0.0f32; n_clusters];
    let mut start = cluster_to_end[0];
    for j in 1..n_clusters {
        let end = cluster_to_end[j];
        let mut p_no_recomb = 1.0f32;
        for k in start..end {
            p_no_recomb *= 1.0f32 - p_recomb.get(k);
        }
        p_clust_recomb[j] = 1.0f32 - p_no_recomb;
        start = end;
    }
    FloatArray::from_floats(&p_clust_recomb)
}

impl MarkerCluster {
    /// `new MarkerCluster(PhaseData phaseData, int sample)`.
    pub fn new(phase_data: &PhaseData, sample: i32) -> Self {
        let sample_phase = phase_data.est_phase().get(sample);
        let cluster_to_end = sample_phase.clust_ends();
        let unph_het_clusters = unph_het_clusters(&sample_phase);
        let p_recomb = p_clust_recomb(&phase_data.p_recomb(), &cluster_to_end);
        MarkerCluster {
            sample_phase,
            cluster_to_end,
            unph_het_clusters,
            p_recomb,
        }
    }

    /// `samplePhase()`.
    pub fn sample_phase(&self) -> &SamplePhase {
        &self.sample_phase
    }

    /// `nClusters()`.
    pub fn n_clusters(&self) -> i32 {
        self.cluster_to_end.len() as i32
    }

    /// `clusterStart(int index)` — inclusive start marker of the cluster.
    pub fn cluster_start(&self, index: i32) -> i32 {
        if index == 0 {
            0
        } else {
            self.cluster_to_end[(index - 1) as usize]
        }
    }

    /// `clusterEnd(int index)` — exclusive end marker of the cluster.
    pub fn cluster_end(&self, index: i32) -> i32 {
        self.cluster_to_end[index as usize]
    }

    /// `pRecomb()` — per-cluster probability of transitioning to a random HMM state.
    pub fn p_recomb(&self) -> &FloatArray {
        &self.p_recomb
    }

    /// `unphasedHetClusters()` — increasing cluster indices containing an unphased het.
    pub fn unphased_het_clusters(&self) -> &WrappedIntArray {
        &self.unph_het_clusters
    }

    /// `isUnphasedHet(int cluster)`.
    pub fn is_unphased_het(&self, cluster: i32) -> bool {
        self.sample_phase.clust_type(cluster) == ClustType::UnphasedHet
    }

    /// `isPhasedHet(int cluster)`.
    pub fn is_phased_het(&self, cluster: i32) -> bool {
        self.sample_phase.clust_type(cluster) == ClustType::PhasedHet
    }

    /// `isHet(int cluster)`.
    pub fn is_het(&self, cluster: i32) -> bool {
        let ct = self.sample_phase.clust_type(cluster);
        ct == ClustType::UnphasedHet || ct == ClustType::PhasedHet
    }

    /// `isMissingGT(int cluster)`.
    pub fn is_missing_gt(&self, cluster: i32) -> bool {
        self.sample_phase.clust_type(cluster) == ClustType::MissingGt
    }

    /// `isMaskedHet(int cluster)`.
    pub fn is_masked_het(&self, cluster: i32) -> bool {
        self.sample_phase.clust_type(cluster) == ClustType::MaskedHet
    }

    /// `isMissingGtOrMaskedHet(int cluster)`.
    pub fn is_missing_gt_or_masked_het(&self, cluster: i32) -> bool {
        let ct = self.sample_phase.clust_type(cluster);
        ct == ClustType::MissingGt || ct == ClustType::MaskedHet
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::IntArray;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, GT, HEADER_PREFIX,
    };
    use std::io::Write;
    use std::rc::Rc;

    use super::super::FixedPhaseData;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_mc_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..3 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 3],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let recs: Vec<Rc<dyn GTRec>> = (0..4)
            .map(|i| {
                let pos = 1_000_000 + i * 1_000_000;
                let line = format!("chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1");
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(
                    &h, &line, &mp,
                ))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(recs);
        let n = targ.n_markers();
        let indices = MarkerIndices::from_counts(0, n, n);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        let w = Window::new(gen_map, 1, true, indices, None, targ);
        let ped = Pedigree::new(w.targ_gt().samples().clone(), None);
        FixedPhaseData::new(&par, &ped, &w, None)
    }

    #[test]
    fn clusters_cover_markers_and_precomb_is_valid() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let n_markers = pd.fpd().stage1_targ_gt().n_markers();

        // sample 1 (0|1 at every marker) has hets -> multiple clusters
        let mc = MarkerCluster::new(&pd, 1);
        let n_clust = mc.n_clusters();
        assert!(n_clust >= 1);
        assert_eq!(mc.cluster_start(0), 0);
        assert_eq!(mc.cluster_end(n_clust - 1), n_markers);
        // contiguous, non-decreasing boundaries
        for c in 1..n_clust {
            assert_eq!(mc.cluster_start(c), mc.cluster_end(c - 1));
        }
        // p_recomb has one entry per cluster, in [0, 1], with index 0 == 0.0
        let pr = mc.p_recomb();
        assert_eq!(pr.size(), n_clust);
        assert_eq!(pr.get(0), 0.0);
        for c in 0..n_clust {
            assert!((0.0..=1.0).contains(&pr.get(c)));
        }

        // sample 0 (0|0 homozygous) has no unphased hets
        let mc0 = MarkerCluster::new(&pd, 0);
        assert_eq!(mc0.unphased_het_clusters().size(), 0);
        for c in 0..mc0.n_clusters() {
            assert!(!mc0.is_het(c));
        }
    }
}
