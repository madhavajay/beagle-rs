//! Port of `blbutil/FileIt.java` (trait) and `blbutil/InputIt.java` (line iterator).
//!
//! `InputIt` yields the lines of a plain, GZIP, or BGZIP text file. For input, only the
//! *line sequence* must match the Java reference (decompressed bytes are unambiguous), so
//! this reads via `flate2::MultiGzDecoder` (which, like Java's `GZIPInputStream`, reads
//! concatenated gzip members = BGZF blocks). The parallel `BGZipIt`/`BlockLineReader`
//! Java classes are a read-throughput optimization with identical output and are
//! subsumed here (see `from_bgzip_file`).

use crate::blbutil::Utilities;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Port of `blbutil/FileIt.java` — an iterator over file elements that knows its file.
pub trait FileIt: Iterator {
    /// The source file, or `None` for stdin / unknown.
    fn file(&self) -> Option<&Path>;
}

/// Port of `blbutil/InputIt.java` — a buffered line iterator.
pub struct InputIt {
    file: Option<PathBuf>,
    reader: Box<dyn BufRead>,
}

fn has_gzip_suffix(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|e| e.to_str()),
        Some("gz") | Some("bgz")
    )
}

impl InputIt {
    /// Builds an `InputIt` from an arbitrary buffered reader (used for stdin / testing).
    pub fn from_reader(reader: Box<dyn BufRead>, file: Option<PathBuf>) -> Self {
        InputIt { file, reader }
    }

    /// `InputIt.fromTextFile(File)`.
    pub fn from_text_file(path: &Path) -> Self {
        let f = Self::open(path);
        InputIt {
            file: Some(path.to_path_buf()),
            reader: Box::new(BufReader::new(f)),
        }
    }

    /// `InputIt.fromGzipFile(File)` — GZIP/BGZIP-decompresses if the name ends in
    /// `.gz`/`.bgz`, else reads plain text.
    pub fn from_gzip_file(path: &Path) -> Self {
        let f = Self::open(path);
        let reader: Box<dyn BufRead> = if has_gzip_suffix(path) {
            Box::new(BufReader::new(MultiGzDecoder::new(f)))
        } else {
            Box::new(BufReader::new(f))
        };
        InputIt {
            file: Some(path.to_path_buf()),
            reader,
        }
    }

    /// `InputIt.fromBGZipFile(File, int)` — same observable line sequence as
    /// `from_gzip_file` (the `nBufferedBlocks` parallel-prefetch hint has no effect on
    /// output).
    pub fn from_bgzip_file(path: &Path, _n_buffered_blocks: i32) -> Self {
        Self::from_gzip_file(path)
    }

    fn open(path: &Path) -> File {
        match File::open(path) {
            Ok(f) => f,
            Err(e) => Utilities::exit(&format!("Error opening {}: {e}", path.display())),
        }
    }
}

impl Iterator for InputIt {
    type Item = String;

    fn next(&mut self) -> Option<String> {
        let mut buf = String::new();
        match self.reader.read_line(&mut buf) {
            Ok(0) => None,
            Ok(_) => {
                // Match BufferedReader.readLine(): strip a trailing \n and a preceding \r.
                if buf.ends_with('\n') {
                    buf.pop();
                    if buf.ends_with('\r') {
                        buf.pop();
                    }
                }
                Some(buf)
            }
            Err(e) => {
                let name = self
                    .file
                    .as_deref()
                    .map_or("<stream>".into(), |p| p.display().to_string());
                Utilities::exit(&format!("Error reading {name}: {e}"))
            }
        }
    }
}

impl FileIt for InputIt {
    fn file(&self) -> Option<&Path> {
        self.file.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blbutil::BgzipOutputStream;
    use std::io::{Cursor, Write};

    fn lines_from(bytes: Vec<u8>) -> Vec<String> {
        InputIt::from_reader(Box::new(BufReader::new(Cursor::new(bytes))), None).collect()
    }

    #[test]
    fn reads_plain_lines_stripping_terminators() {
        let lines = lines_from(b"line1\nline2\r\nline3\nno-newline".to_vec());
        assert_eq!(lines, vec!["line1", "line2", "line3", "no-newline"]);
    }

    #[test]
    fn reads_bgzip_lines() {
        let text = "##header\n#CHROM\tPOS\nchr1\t100\nchr1\t200\n";
        let mut s = BgzipOutputStream::new(Vec::new(), true);
        s.write_all(text.as_bytes()).unwrap();
        let gz = s.close().unwrap();
        let reader = BufReader::new(MultiGzDecoder::new(Cursor::new(gz)));
        let lines: Vec<String> = InputIt::from_reader(Box::new(reader), None).collect();
        assert_eq!(
            lines,
            vec!["##header", "#CHROM\tPOS", "chr1\t100", "chr1\t200"]
        );
    }

    #[test]
    fn reads_real_beagle_bgzip_fixture() {
        // Reads Beagle's own BGZF output (the reference panel) and counts lines.
        let path = format!(
            "{}/../../fixtures/official/inputs/ref.27Feb25.75f.vcf.gz",
            env!("CARGO_MANIFEST_DIR")
        );
        let n = InputIt::from_gzip_file(Path::new(&path)).count();
        assert_eq!(n, 1386); // 30 meta-info lines + 1 header + 1356 marker records (was wc -l)
    }
}
