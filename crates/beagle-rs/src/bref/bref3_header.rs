//! Port of `bref/Bref3Header.java` — the header of a bref3 file: magic number, encoding
//! program name, and sample identifiers, with the sample/haplotype filter mappings.

use std::io::{self, Read};
use std::path::Path;

use crate::blbutil::{consts, Filter, Utilities};
use crate::jdk_io::DataIn;
use crate::vcf::Samples;

use super::{read_string_array, MAGIC_NUMBER_V3};

/// Port of `bref/Bref3Header.java`.
pub struct Bref3Header {
    program: String,
    sample_ids: Vec<String>,
    filtered_hap_indices: Vec<i32>,
    inv_filtered_hap_indices: Vec<i32>,
    samples: Samples,
}

fn check_magic_number(magic_number: i32) {
    if magic_number != MAGIC_NUMBER_V3 {
        Utilities::exit(&format!(
            "ERROR: Unrecognized input file.  Was input file created {}\
             with a different version of the bref program?",
            consts::NL
        ));
    }
}

fn read_program_and_ids<R: Read>(di: &mut DataIn<R>) -> (String, Vec<String>) {
    let result = (|| -> io::Result<(String, Vec<String>)> {
        check_magic_number(di.read_int()?);
        let program = di.read_utf()?;
        let ids = read_string_array(di)?.expect("sample id array");
        Ok((program, ids))
    })();
    result.unwrap_or_else(|e| Utilities::exit(&format!("Error reading file: {e}")))
}

fn filtered_sample_indices(ids: &[String], filter: &Filter<String>) -> Vec<i32> {
    (0..ids.len() as i32)
        .filter(|&j| filter.accept(&ids[j as usize]))
        .collect()
}

fn no_samples_error(source: &str) -> ! {
    Utilities::exit(&format!(
        "{nl}Error      :  All samples in the bref3 file have been excluded\
         {nl}Bref3 file :  {source}",
        nl = consts::NL
    ))
}

/// `hapIndices(int[] sampIndices)` — `[2s, 2s+1]` for each sample index `s`.
fn hap_indices(samp_indices: &[i32]) -> Vec<i32> {
    assert!(samp_indices.len() < (1 << 30), "{}", samp_indices.len());
    let mut hap_indices = Vec::with_capacity(samp_indices.len() << 1);
    for &s in samp_indices {
        let hap1 = s << 1;
        hap_indices.push(hap1);
        hap_indices.push(hap1 | 0b1);
    }
    hap_indices
}

/// `invArray(int[] hapIndices, int size)` — inverse map (value `-1` where absent).
fn inv_array(hap_indices: &[i32], size: i32) -> Vec<i32> {
    assert!(size >= 0, "{size}");
    let mut inverse = vec![-1i32; size as usize];
    for (j, &h) in hap_indices.iter().enumerate() {
        assert!(inverse[h as usize] < 0, "duplicate array value: {h}");
        inverse[h as usize] = j as i32;
    }
    inverse
}

fn build_samples(sample_ids: &[String], sample_indices: &[i32]) -> Samples {
    let ids: Vec<String> = sample_indices
        .iter()
        .map(|&i| sample_ids[i as usize].clone())
        .collect();
    let is_diploid = vec![true; ids.len()];
    Samples::new(&ids, &is_diploid)
}

impl Bref3Header {
    /// `new Bref3Header(File source, DataInput dataIn, Filter<String> sampleFilter)`.
    pub fn new<R: Read>(
        source: Option<&Path>,
        di: &mut DataIn<R>,
        sample_filter: &Filter<String>,
    ) -> Self {
        let (program, sample_ids) = read_program_and_ids(di);
        let filtered_sample_indices = filtered_sample_indices(&sample_ids, sample_filter);
        if filtered_sample_indices.is_empty() {
            let src = source.map_or_else(|| "stdin".to_string(), |p| p.display().to_string());
            no_samples_error(&src);
        }
        let filtered_hap_indices = hap_indices(&filtered_sample_indices);
        let inv_filtered_hap_indices =
            inv_array(&filtered_hap_indices, (sample_ids.len() << 1) as i32);
        let samples = build_samples(&sample_ids, &filtered_sample_indices);
        Bref3Header {
            program,
            sample_ids,
            filtered_hap_indices,
            inv_filtered_hap_indices,
            samples,
        }
    }

    /// `program()`.
    pub fn program(&self) -> &str {
        &self.program
    }

    /// `unfilteredSampleIds()`.
    pub fn unfiltered_sample_ids(&self) -> Vec<String> {
        self.sample_ids.clone()
    }

    /// `filteredHapIndices()`.
    pub fn filtered_hap_indices(&self) -> Vec<i32> {
        self.filtered_hap_indices.clone()
    }

    /// `invfilteredHapIndices()`.
    pub fn inv_filtered_hap_indices(&self) -> Vec<i32> {
        self.inv_filtered_hap_indices.clone()
    }

    /// `samples()`.
    pub fn samples(&self) -> &Samples {
        &self.samples
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jdk_io::DataOut;
    use std::io::Cursor;

    fn header_bytes(program: &str, ids: &[&str]) -> Vec<u8> {
        let mut out = DataOut::new(Vec::new());
        out.write_int(MAGIC_NUMBER_V3).unwrap();
        out.write_utf(program).unwrap();
        out.write_int(ids.len() as i32).unwrap();
        for id in ids {
            out.write_utf(id).unwrap();
        }
        out.into_inner()
    }

    #[test]
    fn reads_header_accept_all() {
        let bytes = header_bytes("beagle.bref3", &["S0", "S1", "S2"]);
        let mut di = DataIn::new(Cursor::new(bytes));
        let h = Bref3Header::new(None, &mut di, &Filter::accept_all());
        assert_eq!(h.program(), "beagle.bref3");
        assert_eq!(h.unfiltered_sample_ids(), vec!["S0", "S1", "S2"]);
        assert_eq!(h.samples().size(), 3);
        assert_eq!(h.filtered_hap_indices(), vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(h.inv_filtered_hap_indices(), vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn applies_sample_filter() {
        let bytes = header_bytes("prog", &["S0", "S1", "S2"]);
        let mut di = DataIn::new(Cursor::new(bytes));
        // include only S1
        let filter = Filter::include(["S1".to_string()]);
        let h = Bref3Header::new(None, &mut di, &filter);
        assert_eq!(h.samples().size(), 1);
        assert_eq!(h.samples().id(0), "S1");
        // S1 is sample index 1 -> haps [2, 3]
        assert_eq!(h.filtered_hap_indices(), vec![2, 3]);
        // inverse over 6 unfiltered haps
        assert_eq!(h.inv_filtered_hap_indices(), vec![-1, -1, 0, 1, -1, -1]);
    }
}
