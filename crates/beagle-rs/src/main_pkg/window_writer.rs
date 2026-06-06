//! Port of `main/WindowWriter.java` — writes the BGZIP-compressed output VCF
//! (`<out>.vcf.gz`) in independently-compressed blocks.
//!
//! Java compresses each step-range of records into its own `ByteArrayOutputStream` via a
//! `BGZIPOutputStream(baos, false)` (no EOF marker), collects the resulting byte arrays in
//! encounter order with a parallel `IntStream` (`toArray()` preserves index order), and
//! appends them to the file in order. The Rust port runs the per-range compression as a
//! sequential, index-ordered loop — byte-identical, since output order is preserved either
//! way. The final EOF marker is written by [`WindowWriter::close`].
//!
//! Wall-clock note: `printPhased` ends with a `System.out.println("WindWriter tot: ...")`
//! debug line (stdout only — never the `.log` or `.vcf.gz`); it is reproduced faithfully and
//! is a documented stdout normalization line.

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::time::Instant;

use crate::blbutil::{BgzipOutputStream, Utilities};
use crate::imp::{ImpData, ImputedVcfWriter, StateProbs};
use crate::ints::UnsignedByteArray;
use crate::phase::Stage2Haps;
use crate::vcf::{append_records_gt, append_records_recs, write_meta_lines, GTRec, Samples, GT};

use super::{Par, PROGRAM};

/// Port of `main/WindowWriter.java`.
pub struct WindowWriter {
    samples: Samples,
    out_prefix: String,
    vcf_out_file: PathBuf,
}

impl WindowWriter {
    /// `WindowWriter(Par par, Samples samples)` — creates `<out>.vcf.gz` and writes the
    /// (BGZIP-compressed) VCF meta-information and header lines.
    pub fn new(par: &Par, samples: Samples) -> Self {
        let out_prefix = par.out().to_string();
        let vcf_out_file = PathBuf::from(format!("{out_prefix}.vcf.gz"));

        let ds = true;
        let ap = par.ap();
        let gp = par.gp();
        let gl = false;
        let mut bgzip = BgzipOutputStream::new(Vec::new(), false);
        write_meta_lines(&samples.ids(), Some(PROGRAM), ds, ap, gp, gl, &mut bgzip)
            .expect("write VCF meta lines");
        let bytes = bgzip.close().expect("close BGZIP meta block");

        match File::create(&vcf_out_file) {
            Ok(mut fos) => {
                if let Err(e) = fos.write_all(&bytes) {
                    Self::file_output_error(&vcf_out_file, &e);
                }
            }
            Err(e) => Self::file_output_error(&vcf_out_file, &e),
        }

        WindowWriter {
            samples,
            out_prefix,
            vcf_out_file,
        }
    }

    /// `outPrefix()`.
    pub fn out_prefix(&self) -> &str {
        &self.out_prefix
    }

    /// `samples()`.
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// `printImputed(ImpData impData, int start, int end, AtomicReferenceArray<StateProbs>)`.
    pub fn print_imputed(
        &mut self,
        imp_data: &ImpData,
        start: i32,
        end: i32,
        state_probs: &[Box<dyn StateProbs>],
    ) {
        Self::check_interval(start, end, imp_data.ref_gt().n_markers());
        assert!(
            state_probs.len() as i32 == imp_data.n_targ_haps(),
            "inconsistent data:"
        );
        let output: Vec<UnsignedByteArray> = (0..imp_data.n_clusters())
            .map(|c| Self::imputed_to_byte_array(imp_data, start, end, state_probs, c))
            .collect();
        Self::append(&output, &self.vcf_out_file);
    }

    fn imputed_to_byte_array(
        imp_data: &ImpData,
        ref_start: i32,
        ref_end: i32,
        state_probs: &[Box<dyn StateProbs>],
        m: i32,
    ) -> UnsignedByteArray {
        let mut ivw = ImputedVcfWriter::new(imp_data, ref_start, ref_end, m);
        let mut bgzip = BgzipOutputStream::new(Vec::new(), false);
        ivw.append_records(state_probs, &mut bgzip);
        UnsignedByteArray::from_baos(bgzip.close().expect("close BGZIP block"))
    }

    /// `printPhased(GT phasedTarg, int start, int end)`.
    pub fn print_phased_gt(&mut self, phased_targ: &dyn GT, start: i32, end: i32) {
        Self::check_interval(start, end, phased_targ.n_markers());
        let t0 = Instant::now();
        let block_size = 50000;
        let step_size = 100;
        let block_ends = Self::ends(start, end, block_size);
        for j in 1..block_ends.len() {
            let step_ends = Self::ends(block_ends[j - 1], block_ends[j], step_size);
            let output: Vec<UnsignedByteArray> = (1..step_ends.len())
                .map(|i| Self::gt_to_byte_array(phased_targ, step_ends[i - 1], step_ends[i]))
                .collect();
            Self::append(&output, &self.vcf_out_file);
        }
        let t1 = t0.elapsed();
        println!(
            "WindWriter tot: {}",
            Utilities::elapsed_nanos(t1.as_nanos() as i64)
        );
    }

    /// `printPhased(Stage2Haps stage2Haps, int start, int end)`.
    pub fn print_phased_stage2(&mut self, stage2_haps: &Stage2Haps, start: i32, end: i32) {
        let targ_gt = stage2_haps.fpd().targ_gt();
        Self::check_interval(start, end, targ_gt.n_markers());
        let t0 = Instant::now();
        let block_size = 20000;
        let step_size = 50;
        let block_ends = Self::ends(start, end, block_size);
        for j in 1..block_ends.len() {
            let block_start = block_ends[j - 1];
            let block_end = block_ends[j];
            let rec_block = stage2_haps.to_gt_recs(block_start, block_end);
            let step_ends = Self::ends(block_start, block_end, step_size);
            let output: Vec<UnsignedByteArray> = (1..step_ends.len())
                .map(|i| {
                    let from = (step_ends[i - 1] - block_start) as usize;
                    let to = (step_ends[i] - block_start) as usize;
                    Self::recs_to_byte_array(&self.samples, &rec_block[from..to])
                })
                .collect();
            Self::append(&output, &self.vcf_out_file);
        }
        let t1 = t0.elapsed();
        println!(
            "WindWriter tot: {}",
            Utilities::elapsed_nanos(t1.as_nanos() as i64)
        );
    }

    fn check_interval(start: i32, end: i32, n_markers: i32) {
        assert!(
            !(start < 0 || end > n_markers || end < start),
            "start={start} end={end} nMarkers={n_markers}"
        );
    }

    /// `ends(int start, int end, int step)` — block/step boundaries `[start, start+step, …, end]`.
    fn ends(start: i32, end: i32, step: i32) -> Vec<i32> {
        let mut starts = Vec::with_capacity((2 + (end - start) / step) as usize);
        let mut m = start;
        while m < end {
            starts.push(m);
            m += step;
        }
        starts.push(end);
        starts
    }

    fn gt_to_byte_array(phased_targ: &dyn GT, start: i32, end: i32) -> UnsignedByteArray {
        let mut bgzip = BgzipOutputStream::new(Vec::new(), false);
        append_records_gt(phased_targ, start, end, &mut bgzip).expect("append GT records");
        UnsignedByteArray::from_baos(bgzip.close().expect("close BGZIP block"))
    }

    fn recs_to_byte_array(samples: &Samples, phased_recs: &[Rc<dyn GTRec>]) -> UnsignedByteArray {
        let refs: Vec<&dyn GTRec> = phased_recs.iter().map(|r| r.as_ref()).collect();
        let mut bgzip = BgzipOutputStream::new(Vec::new(), false);
        append_records_recs(samples, &refs, &mut bgzip).expect("append rec records");
        UnsignedByteArray::from_baos(bgzip.close().expect("close BGZIP block"))
    }

    fn append(output: &[UnsignedByteArray], out_file: &Path) {
        match OpenOptions::new().append(true).open(out_file) {
            Ok(f) => {
                let mut bos = BufWriter::new(f);
                for uba in output {
                    if let Err(e) = bos.write_all(uba.as_bytes()) {
                        Self::file_output_error(out_file, &e);
                    }
                }
                if let Err(e) = bos.flush() {
                    Self::file_output_error(out_file, &e);
                }
            }
            Err(e) => Self::file_output_error(out_file, &e),
        }
    }

    /// `close()` — appends the empty BGZIP block (the BGZF EOF marker) to `<out>.vcf.gz`.
    pub fn close(&mut self) {
        match OpenOptions::new().append(true).open(&self.vcf_out_file) {
            Ok(f) => {
                let bgzip = BgzipOutputStream::new(BufWriter::new(f), true);
                if let Err(e) = bgzip.close() {
                    Utilities::exit(&format!(
                        "Error closing file: {} : {e}",
                        self.vcf_out_file.display()
                    ));
                }
            }
            Err(e) => Utilities::exit(&format!(
                "Error closing file: {} : {e}",
                self.vcf_out_file.display()
            )),
        }
    }

    fn file_output_error(file: &Path, e: &std::io::Error) -> ! {
        Utilities::exit(&format!("Error writing to file: {} : {e}", file.display()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ends_partitions_interval() {
        // start == end -> single element, so no step ranges are produced
        assert_eq!(WindowWriter::ends(5, 5, 100), vec![5]);
        // exact multiples
        assert_eq!(WindowWriter::ends(0, 250, 100), vec![0, 100, 200, 250]);
        // step larger than the interval
        assert_eq!(WindowWriter::ends(0, 250, 50000), vec![0, 250]);
        // boundary lands exactly on `end`
        assert_eq!(WindowWriter::ends(0, 200, 100), vec![0, 100, 200]);
    }
}
