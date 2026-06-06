//! Port of `imp/ImpData.java` — the immutable per-window input data for imputing ungenotyped
//! markers: target marker clusters, cluster genetic-map positions, per-cluster error and
//! recombination probabilities, reference-marker interpolation weights, and the coded
//! reference+target allele sequences per cluster.

use std::rc::Rc;

use crate::ints::IndexArray;
use crate::main_pkg::Par;
use crate::vcf::{GeneticMap, MarkerIndices, Markers, RefGT, Window, GT};

use super::HaplotypeCoder;

const MIN_CM_DIST: f64 = 1e-7;

/// Port of `imp/ImpData.java`.
pub struct ImpData {
    par: Par,
    marker_indices: MarkerIndices,
    ref_gt: RefGT,
    phased_targ: Rc<dyn GT>,
    targ_clust_start_end: Vec<i32>,
    ref_cluster_start: Vec<i32>,
    ref_cluster_end: Vec<i32>,
    hap_to_seq: Vec<IndexArray>,
    err_prob: Vec<f32>,
    pos: Vec<f64>,
    p_recomb: Vec<f32>,
    weight: Vec<f32>,
    n_clusters: i32,
    n_ref_haps: i32,
    n_targ_haps: i32,
    n_input_targ_haps: i32,
    n_haps: i32,
}

fn cum_pos(markers: &Markers, map: &dyn GeneticMap) -> Vec<f64> {
    let n = markers.size() as usize;
    let mut cum_pos = vec![0.0f64; n];
    let mut last_gen_pos = map.gen_pos_marker(markers.marker(0));
    for j in 1..n {
        let gen_pos = map.gen_pos_marker(markers.marker(j as i32));
        let gen_dist = (gen_pos - last_gen_pos).abs().max(MIN_CM_DIST);
        cum_pos[j] = cum_pos[j - 1] + gen_dist;
        last_gen_pos = gen_pos;
    }
    cum_pos
}

fn targ_block_end(ref_gt: &RefGT, targ_to_ref: &[i32]) -> Vec<i32> {
    let mut int_list: Vec<i32> = Vec::with_capacity(targ_to_ref.len() / 4);
    let mut last_key: Option<usize> = None; // None == Java's null `lastHap2Seq`
    for (j, &ref_index) in targ_to_ref.iter().enumerate() {
        let rec = ref_gt.get(ref_index);
        if !rec.is_allele_coded() {
            let key = rec.seq_block_key(); // Some(addr) for seq-coded
            if key != last_key {
                if last_key.is_some() {
                    int_list.push(j as i32);
                }
                last_key = key;
            }
        }
    }
    int_list.push(targ_to_ref.len() as i32);
    int_list
}

fn targ_clust_start_end(raw_pos: &[f64], targ_block_end: &[i32], cluster_dist: f32) -> Vec<i32> {
    let mut clust_start_end = vec![0i32; raw_pos.len() + 1];
    let mut size = 1usize; // clust_start_end[0] = 0
    for &block_end in targ_block_end {
        let clust_start = clust_start_end[size - 1];
        let mut start_pos = raw_pos[clust_start as usize];
        for m in (clust_start + 1)..block_end {
            let pos = raw_pos[m as usize];
            if (pos - start_pos) > cluster_dist as f64 {
                clust_start_end[size] = m;
                size += 1;
                start_pos = pos;
            }
        }
        clust_start_end[size] = block_end;
        size += 1;
        // Java also writes targBlockEnd[j] = size-2 here, but that value is never read.
    }
    clust_start_end.truncate(size);
    clust_start_end
}

fn mid_pos(pos: &[f64], start_end: &[i32]) -> Vec<f64> {
    (1..start_end.len())
        .map(|j| (pos[start_end[j - 1] as usize] + pos[(start_end[j] - 1) as usize]) / 2.0)
        .collect()
}

fn build_hap_to_seq(
    restrict_ref: RefGT,
    phased_targ: Rc<dyn GT>,
    targ_start_end: &[i32],
) -> Vec<IndexArray> {
    let coder = HaplotypeCoder::new(restrict_ref, phased_targ);
    (1..targ_start_end.len())
        .map(|j| coder.run(targ_start_end[j - 1], targ_start_end[j]))
        .collect()
}

fn err(err_rate: f32, start_end: &[i32]) -> Vec<f32> {
    let max_err_prob = 0.5f32;
    let mut err = vec![0f32; start_end.len() - 1];
    for j in 0..err.len() {
        err[j] = err_rate * (start_end[j + 1] - start_end[j]) as f32;
        if err[j] > max_err_prob {
            err[j] = max_err_prob;
        }
    }
    err
}

fn ref_clust_start(clust_start_end: &[i32], targ_to_ref: &[i32]) -> Vec<i32> {
    (0..clust_start_end.len() - 1)
        .map(|j| targ_to_ref[clust_start_end[j] as usize])
        .collect()
}

fn ref_clust_end(clust_start_end: &[i32], targ_to_ref: &[i32]) -> Vec<i32> {
    (1..clust_start_end.len())
        .map(|j| targ_to_ref[(clust_start_end[j] - 1) as usize] + 1)
        .collect()
}

fn p_recomb_arr(ne: f32, n_haps: i32, pos: &[f64]) -> Vec<f32> {
    let mut p_recomb = vec![0f32; pos.len()];
    let c = -(0.04f64 * ne as f64 / n_haps as f64); // 0.04 = 4/(100 cM/M)
    for j in 1..p_recomb.len() {
        p_recomb[j] = -((c * (pos[j] - pos[j - 1])).exp_m1()) as f32;
    }
    p_recomb
}

fn wts(
    ref_markers: &Markers,
    ref_cluster_start: &[i32],
    ref_cluster_end: &[i32],
    map: &dyn GeneticMap,
) -> Vec<f32> {
    let cum_pos = cum_pos(ref_markers, map);
    let n_targ_markers_m1 = ref_cluster_start.len() - 1;
    let mut wts = vec![0f32; cum_pos.len()];
    wts[0..ref_cluster_start[0] as usize].fill(f32::NAN);
    for j in 0..n_targ_markers_m1 {
        let start = ref_cluster_start[j];
        let end = ref_cluster_end[j];
        let next_start = ref_cluster_start[j + 1];
        let total_length = cum_pos[next_start as usize] - cum_pos[(end - 1) as usize];
        wts[start as usize..end as usize].fill(f32::NAN);
        for m in end..next_start {
            wts[m as usize] =
                ((cum_pos[next_start as usize] - cum_pos[m as usize]) / total_length) as f32;
        }
    }
    wts[ref_cluster_start[n_targ_markers_m1] as usize..ref_markers.size() as usize].fill(f32::NAN);
    wts
}

impl ImpData {
    /// `new ImpData(Par par, Window window, GT phasedTarg, GeneticMap map)`.
    pub fn new(par: &Par, window: &Window, phased_targ: Rc<dyn GT>, map: &dyn GeneticMap) -> Self {
        assert!(
            window.targ_gt().markers() == phased_targ.markers(),
            "inconsistent markers"
        );
        assert!(
            window.targ_gt().samples() == phased_targ.samples(),
            "inconsistent samples"
        );
        assert!(phased_targ.is_phased(), "unphased data");
        let marker_indices = window.indices().clone();
        let ref_gt = window.ref_gt().expect("ref present").clone();
        let targ_to_ref = marker_indices.targ_marker_to_marker();
        let targ_pos = cum_pos(phased_targ.markers(), map);
        let block_end = targ_block_end(&ref_gt, &targ_to_ref);
        let targ_clust_start_end = targ_clust_start_end(&targ_pos, &block_end, par.cluster());
        let pos = mid_pos(&targ_pos, &targ_clust_start_end);
        let n_clusters = targ_clust_start_end.len() as i32 - 1;
        let n_ref_haps = ref_gt.n_haps();
        let n_targ_haps = phased_targ.n_haps();
        let targ_samples = phased_targ.samples();
        let n_input_targ_haps: i32 = (0..targ_samples.size())
            .map(|j| if targ_samples.is_diploid(j) { 2 } else { 1 })
            .sum();
        let n_haps = n_ref_haps + n_targ_haps;
        let hap_to_seq = build_hap_to_seq(
            window
                .restrict_ref_gt()
                .expect("restrict ref present")
                .clone(),
            phased_targ.clone(),
            &targ_clust_start_end,
        );
        let ref_cluster_start = ref_clust_start(&targ_clust_start_end, &targ_to_ref);
        let ref_cluster_end = ref_clust_end(&targ_clust_start_end, &targ_to_ref);
        let err_prob = err(par.err(n_haps), &targ_clust_start_end);
        let p_recomb = p_recomb_arr(par.ne(), ref_gt.n_haps(), &pos);
        let weight = wts(ref_gt.markers(), &ref_cluster_start, &ref_cluster_end, map);
        ImpData {
            par: par.clone(),
            marker_indices,
            ref_gt,
            phased_targ,
            targ_clust_start_end,
            ref_cluster_start,
            ref_cluster_end,
            hap_to_seq,
            err_prob,
            pos,
            p_recomb,
            weight,
            n_clusters,
            n_ref_haps,
            n_targ_haps,
            n_input_targ_haps,
            n_haps,
        }
    }

    /// `par()`.
    pub fn par(&self) -> &Par {
        &self.par
    }
    /// `markerIndices()`.
    pub fn marker_indices(&self) -> &MarkerIndices {
        &self.marker_indices
    }
    /// `refGT()`.
    pub fn ref_gt(&self) -> &RefGT {
        &self.ref_gt
    }
    /// `targGT()`.
    pub fn targ_gt(&self) -> &Rc<dyn GT> {
        &self.phased_targ
    }
    /// `targClusterStart(int cluster)`.
    pub fn targ_cluster_start(&self, cluster: i32) -> i32 {
        assert!(cluster < self.n_clusters, "{cluster}");
        self.targ_clust_start_end[cluster as usize]
    }
    /// `targClusterEnd(int cluster)`.
    pub fn targ_cluster_end(&self, cluster: i32) -> i32 {
        assert!(cluster >= 0, "{cluster}");
        self.targ_clust_start_end[(cluster + 1) as usize]
    }
    /// `refClusterStart(int cluster)`.
    pub fn ref_cluster_start(&self, cluster: i32) -> i32 {
        self.ref_cluster_start[cluster as usize]
    }
    /// `refClusterEnd(int cluster)`.
    pub fn ref_cluster_end(&self, cluster: i32) -> i32 {
        self.ref_cluster_end[cluster as usize]
    }
    /// `nClusters()`.
    pub fn n_clusters(&self) -> i32 {
        self.n_clusters
    }
    /// `targSamples()` / `nTargSamples()` helpers.
    pub fn n_targ_samples(&self) -> i32 {
        self.phased_targ.n_samples()
    }
    /// `nHaps()`.
    pub fn n_haps(&self) -> i32 {
        self.n_haps
    }
    /// `nRefHaps()`.
    pub fn n_ref_haps(&self) -> i32 {
        self.n_ref_haps
    }
    /// `nTargHaps()`.
    pub fn n_targ_haps(&self) -> i32 {
        self.n_targ_haps
    }
    /// `nInputTargHaps()`.
    pub fn n_input_targ_haps(&self) -> i32 {
        self.n_input_targ_haps
    }
    /// `allele(int cluster, int hap)`.
    pub fn allele(&self, cluster: i32, hap: i32) -> i32 {
        self.hap_to_seq[cluster as usize].int_array().get(hap)
    }
    /// `hapToSeq(int cluster)`.
    pub fn hap_to_seq(&self, cluster: i32) -> &IndexArray {
        &self.hap_to_seq[cluster as usize]
    }
    /// `errProb(int cluster)`.
    pub fn err_prob(&self, cluster: i32) -> f32 {
        self.err_prob[cluster as usize]
    }
    /// `pos(int cluster)`.
    pub fn pos(&self, cluster: i32) -> f64 {
        self.pos[cluster as usize]
    }
    /// `pos()`.
    pub fn pos_all(&self) -> Vec<f64> {
        self.pos.clone()
    }
    /// `pRecomb(int cluster)`.
    pub fn p_recomb(&self, cluster: i32) -> f32 {
        self.p_recomb[cluster as usize]
    }
    /// `weight(int refMarker)`.
    pub fn weight(&self, ref_marker: i32) -> f64 {
        self.weight[ref_marker as usize] as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{
        allele_ref_gt_rec_from_parser, BasicGT, BasicGTRec, GTRec, MarkerIndices, MarkerParser,
        PositionMap, RefGTRec, VcfHeader, VcfRecGTParser, Window, HEADER_PREFIX,
    };

    fn header(n: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n {
            hdr.push_str(&format!("\tS{s}"));
        }
        VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &vec![true; n],
        )
    }

    #[test]
    fn builds_clusters_and_interpolation_arrays() {
        // ref panel: 4 markers, 3 samples; target: markers 0 and 3 genotyped (1 and 2 imputed)
        let gt = std::env::temp_dir().join("beagle_rs_impdata_gt.vcf");
        let rf = std::env::temp_dir().join("beagle_rs_impdata_ref.vcf");
        std::fs::write(&gt, b"x").unwrap();
        std::fs::write(&rf, b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            format!("ref={}", rf.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
        ]);
        let hr = header(3);
        let ht = header(2);
        let mp = MarkerParser::new(true, true, true, true);
        let ref_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|1\t0|1\t1|0",
        ];
        let ref_recs: Vec<Rc<dyn RefGTRec>> = ref_lines
            .iter()
            .map(|l| {
                Rc::from(allele_ref_gt_rec_from_parser(&VcfRecGTParser::new(
                    &hr, l, &mp,
                )))
            })
            .collect();
        let ref_gt = RefGT::from_recs(ref_recs);

        // target has only markers 100 and 400 (ref markers 0 and 3)
        let targ_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|1",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|0\t1|0",
        ];
        let targ_recs: Vec<Rc<dyn GTRec>> = targ_lines
            .iter()
            .map(|l| {
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(&ht, l, &mp))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(targ_recs);
        // marker indices: 4 total markers, target markers map to ref 0 and 3
        let indices = MarkerIndices::from_in_targ(&[true, false, false, true], 0, 4);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        let window = Window::new(gen_map.clone(), 1, true, indices, Some(ref_gt), targ);

        let phased: Rc<dyn GT> = Rc::new(window.targ_gt().clone());
        let imp = ImpData::new(&par, &window, phased, gen_map.as_ref());

        assert_eq!(imp.n_ref_haps(), 6);
        assert_eq!(imp.n_targ_haps(), 4);
        assert_eq!(imp.n_haps(), 10);
        assert!(imp.n_clusters() >= 1);
        // each cluster's coded alleles cover all haps
        for c in 0..imp.n_clusters() {
            assert!(imp.targ_cluster_start(c) < imp.targ_cluster_end(c));
            for hap in 0..imp.n_haps() {
                let _ = imp.allele(c, hap);
            }
            assert!(imp.err_prob(c) >= 0.0 && imp.err_prob(c) <= 0.5);
            assert!(imp.p_recomb(c) >= 0.0);
        }
    }
}
