//! Port of the Java `bref` package — the bref3 (binary reference format version 3) reader
//! and writers. Built on [`crate::jdk_io`] (`DataInput`/`DataOutput`) for byte-exact output.

use crate::jdk_io::DataIn;
use std::io::{self, Read};

mod as_is_bref3_writer;
mod bref3_header;
mod bref3_it;
mod bref3_reader;
mod bref_block;
mod bref_writer;
mod compress_bref3_writer;
mod seq_coder3;

pub use as_is_bref3_writer::AsIsBref3Writer;
pub use bref3_header::Bref3Header;
pub use bref3_it::Bref3It;
pub use bref3_reader::Bref3Reader;
pub use bref_block::BrefBlock;
pub use bref_writer::BrefWriter;
pub use compress_bref3_writer::CompressBref3Writer;
pub use seq_coder3::{default_max_n_seq, SeqCoder3, COMPRESS_FREQ_THRESHOLD, MAX_NALLELES};

/// `AsIsBref3Writer.END_OF_DATA`.
pub const END_OF_DATA: i32 = 0;
/// `AsIsBref3Writer.END_OF_INDEX`.
pub const END_OF_INDEX: i64 = -999_999_999_999_999;
/// `AsIsBref3Writer.MAGIC_NUMBER_V3`.
pub const MAGIC_NUMBER_V3: i32 = 2_055_763_188;
/// `AsIsBref3Writer.SEQ_CODED`.
pub const SEQ_CODED: i8 = 0;
/// `AsIsBref3Writer.ALLELE_CODED`.
pub const ALLELE_CODED: i8 = 1;

/// `Bref3Reader.readStringArray(DataInput)` — an `int` length then that many modified-UTF-8
/// strings. Returns `None` if the length is negative (Java `null`).
pub(crate) fn read_string_array<R: Read>(di: &mut DataIn<R>) -> io::Result<Option<Vec<String>>> {
    let length = di.read_int()?;
    read_string_array_len(di, length)
}

fn read_string_array_len<R: Read>(
    di: &mut DataIn<R>,
    length: i32,
) -> io::Result<Option<Vec<String>>> {
    if length < 0 {
        Ok(None)
    } else if length == 0 {
        Ok(Some(Vec::new()))
    } else {
        let mut sa = Vec::with_capacity(length as usize);
        for _ in 0..length {
            sa.push(di.read_utf()?);
        }
        Ok(Some(sa))
    }
}
