//! Port of `blbutil/FileUtil.java` (writer side) — buffered file/stdout writers and the
//! BGZIP-compressing writer used for `.vcf.gz` output. Readers are provided by `InputIt`.
//!
//! Java's `PrintWriter` is modelled as `Box<dyn Write>` (always `BufWriter`-buffered). The
//! caller flushes by dropping the writer (Rust `BufWriter` flushes on drop, matching
//! `PrintWriter.close()`); the BGZIP writer additionally needs its `close()` to emit the
//! empty-block EOF marker.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use super::{BgzipOutputStream, Utilities};

/// Port of `blbutil/FileUtil.java` (static methods).
pub struct FileUtil;

impl FileUtil {
    /// `bufferedOutputStream(File file)`.
    pub fn buffered_output_stream(file: &Path) -> Box<dyn Write> {
        match File::create(file) {
            Ok(f) => Box::new(BufWriter::new(f)),
            Err(e) => Utilities::exit(&format!("Error: file not found [{}]: {e}", file.display())),
        }
    }

    /// `printWriter(File file)` — a buffered writer that overwrites `file`.
    pub fn print_writer(file: &Path) -> Box<dyn Write> {
        match File::create(file) {
            Ok(f) => Box::new(BufWriter::new(f)),
            Err(e) => Utilities::exit(&format!("Error opening {}: {e}", file.display())),
        }
    }

    /// `stdOutPrintWriter()`.
    pub fn std_out_print_writer() -> Box<dyn Write> {
        Box::new(BufWriter::new(std::io::stdout()))
    }

    /// `bgzipPrintWriter(File file)` — a BGZIP-compressing writer (its `close()` writes the
    /// empty-block EOF marker).
    pub fn bgzip_print_writer(file: &Path) -> BgzipOutputStream<Box<dyn Write>> {
        BgzipOutputStream::new(Self::buffered_output_stream(file), true)
    }
}
