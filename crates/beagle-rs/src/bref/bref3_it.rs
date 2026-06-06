//! Port of `bref/Bref3It.java` — a `SampleFileIt` whose items are `RefGTRec`s read from a
//! bref3 file (or stdin). Records are produced one block at a time into a FIFO buffer.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use crate::blbutil::{FileIt, Filter, SampleFileIt, Utilities};
use crate::jdk_io::DataIn;
use crate::vcf::{Marker, RefGTRec, Samples};

use super::Bref3Reader;

/// Port of `bref/Bref3It.java`.
pub struct Bref3It {
    bref_file: Option<PathBuf>,
    data_in: DataIn<Box<dyn Read>>,
    bref3_reader: Bref3Reader,
    buffer: VecDeque<Box<dyn RefGTRec>>,
}

impl Bref3It {
    /// `new Bref3It(File brefFile)` — accept-all filters.
    pub fn new(bref_file: Option<&Path>) -> Self {
        Bref3It::with_filters(bref_file, Filter::accept_all(), Filter::accept_all())
    }

    /// `new Bref3It(File brefFile, Filter<String> sampleFilter, Filter<Marker> markerFilter)`.
    pub fn with_filters(
        bref_file: Option<&Path>,
        sample_filter: Filter<String>,
        marker_filter: Filter<Marker>,
    ) -> Self {
        let reader: Box<dyn Read> = match bref_file {
            None => Box::new(BufReader::new(std::io::stdin())),
            Some(f) => match File::open(f) {
                Ok(file) => Box::new(BufReader::new(file)),
                Err(e) => Utilities::exit(&format!("Error opening {}: {e}", f.display())),
            },
        };
        let mut data_in = DataIn::new(reader);
        let bref3_reader =
            Bref3Reader::with_filters(bref_file, &mut data_in, sample_filter, marker_filter);
        let mut bref3_it = Bref3It {
            bref_file: bref_file.map(Path::to_path_buf),
            data_in,
            bref3_reader,
            buffer: VecDeque::with_capacity(500),
        };
        bref3_it
            .bref3_reader
            .read_block(&mut bref3_it.data_in, &mut bref3_it.buffer);
        bref3_it
    }
}

impl Iterator for Bref3It {
    type Item = Box<dyn RefGTRec>;

    fn next(&mut self) -> Option<Box<dyn RefGTRec>> {
        let rec = self.buffer.pop_front()?;
        if self.buffer.is_empty() {
            self.bref3_reader
                .read_block(&mut self.data_in, &mut self.buffer);
        }
        Some(rec)
    }
}

impl FileIt for Bref3It {
    fn file(&self) -> Option<&Path> {
        self.bref_file.as_deref()
    }
}

impl SampleFileIt for Bref3It {
    fn samples(&self) -> &Samples {
        self.bref3_reader.samples()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bref::MAGIC_NUMBER_V3;
    use crate::jdk_io::DataOut;
    use std::io::Write;

    /// Writes a minimal bref3 file: one diploid sample, one seq-coded biallelic SNV.
    fn write_min_bref3() -> PathBuf {
        let mut out = DataOut::new(Vec::new());
        out.write_int(MAGIC_NUMBER_V3).unwrap();
        out.write_utf("test.bref3").unwrap();
        out.write_int(1).unwrap(); // 1 sample id
        out.write_utf("S0").unwrap();
        // --- block ---
        out.write_int(1).unwrap(); // nRecs
        out.write_utf("chr1").unwrap();
        out.write_short(2).unwrap(); // nSeq = 2
                                     // hap->seq map (2 haps, 2 BE bytes each): hap0->seq0, hap1->seq1
        out.write_fully(&[0x00, 0x00, 0x00, 0x01]).unwrap();
        // marker: pos=100, id length 0 (missing), allele_code 1 -> biallelic A/C
        out.write_int(100).unwrap();
        out.write_byte(0).unwrap();
        out.write_byte(1).unwrap();
        // flag 0 = seq-coded; seq->allele map (nSeq=2 bytes): seq0->0, seq1->1
        out.write_byte(0).unwrap();
        out.write_fully(&[0x00, 0x01]).unwrap();
        // end-of-data sentinel
        out.write_int(0).unwrap();

        let path = std::env::temp_dir().join("beagle_rs_min.bref3");
        File::create(&path)
            .unwrap()
            .write_all(&out.into_inner())
            .unwrap();
        path
    }

    #[test]
    fn reads_minimal_bref3_file() {
        let path = write_min_bref3();
        let it = Bref3It::new(Some(&path));
        assert_eq!(it.samples().size(), 1);
        let recs: Vec<Box<dyn RefGTRec>> = it.collect();
        assert_eq!(recs.len(), 1);
        let rec = &recs[0];
        assert_eq!(rec.marker().pos(), 100);
        assert_eq!(rec.size(), 2);
        assert_eq!(rec.get(0), 0); // hap0 -> seq0 -> allele 0 (A)
        assert_eq!(rec.get(1), 1); // hap1 -> seq1 -> allele 1 (C)
        assert!(rec.is_phased());
    }
}
