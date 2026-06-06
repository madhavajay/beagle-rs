//! Port of `blbutil/BGZIPOutputStream.java` — writes the BGZF (blocked gzip) format
//! Beagle uses for `.vcf.gz` output.
//!
//! Byte-for-byte parity with the Java reference requires the raw DEFLATE body to match
//! Java's `Deflater(DEFAULT_COMPRESSION=-1 → level 6, nowrap=true)`, i.e. a zlib-backed
//! deflate (see the crate's `flate2` zlib backend). The block framing (18-byte BGZF
//! header, CRC32, ISIZE), the 65505-byte input block size, the stored-block fallback,
//! and the 28-byte empty EOF block are reproduced exactly.

use flate2::write::DeflateEncoder;
use flate2::{Compression, Crc};
use std::io::{self, Write};

/// `MAX_INPUT_BYTES` = (1<<16) - 31 = 65505: the uncompressed bytes per BGZF block.
pub const MAX_INPUT_BYTES: usize = (1 << 16) - 31;
const NOCOMPRESS_XTRA_BYTES: usize = 5;
const MAX_OUTPUT_BYTES: usize = MAX_INPUT_BYTES + NOCOMPRESS_XTRA_BYTES;

/// Raw DEFLATE (level 6, no zlib wrapper) — matches Java `Deflater(-1, nowrap=true)`.
fn deflate_raw_level6(input: &[u8]) -> Vec<u8> {
    let mut e = DeflateEncoder::new(
        Vec::with_capacity(input.len() / 2 + 128),
        Compression::new(6),
    );
    e.write_all(input).expect("in-memory deflate cannot fail");
    e.finish().expect("in-memory deflate cannot fail")
}

/// `setOutputNoCompression()` — a single stored (BTYPE=00, BFINAL=1) DEFLATE block.
fn stored_block(input: &[u8]) -> Vec<u8> {
    let isize = input.len();
    let lo = (isize & 0xff) as u8;
    let hi = ((isize >> 8) & 0xff) as u8;
    let mut out = Vec::with_capacity(isize + NOCOMPRESS_XTRA_BYTES);
    out.push(1);
    out.push(lo);
    out.push(hi);
    out.push(!lo);
    out.push(!hi);
    out.extend_from_slice(input);
    out
}

fn write_bgzip_header(n_compressed: usize, os: &mut impl Write) -> io::Result<()> {
    let bsize = n_compressed + 25;
    assert!((bsize >> 16) == 0, "{}", n_compressed);
    let header: [u8; 18] = [
        31,
        139, // GZIP magic
        8,   // CM = Deflater.DEFLATED
        4,   // FLG = FEXTRA
        0,
        0,
        0,
        0,   // MTIME = 0
        0,   // XFL
        255, // OS = unknown
        6,
        0, // XLEN = 6
        66,
        67, // BGZF subfield magic "BC"
        2,
        0, // SLEN = 2
        (bsize & 0xff) as u8,
        ((bsize >> 8) & 0xff) as u8, // BSIZE = block size - 1
    ];
    os.write_all(&header)
}

fn write_u32_le(v: u32, os: &mut impl Write) -> io::Result<()> {
    os.write_all(&v.to_le_bytes())
}

fn write_bgzip_block(isize: usize, crc32: u32, body: &[u8], os: &mut impl Write) -> io::Result<()> {
    assert!(isize <= (1 << 16), "{}", isize);
    write_bgzip_header(body.len(), os)?;
    os.write_all(body)?;
    write_u32_le(crc32, os)?;
    write_u32_le(isize as u32, os)
}

/// Compresses one input block (`<= MAX_INPUT_BYTES`) and writes its framed BGZF block.
fn compress_and_write_block(input: &[u8], os: &mut impl Write) -> io::Result<()> {
    let mut crc = Crc::new();
    crc.update(input);
    let crc32 = crc.sum();
    let mut body = deflate_raw_level6(input);
    if body.len() > MAX_OUTPUT_BYTES {
        body = stored_block(input);
    }
    write_bgzip_block(input.len(), crc32, &body, os)
}

/// Writes a single empty BGZF block (the 28-byte BGZF EOF marker).
pub fn write_empty_block(os: &mut impl Write) -> io::Result<()> {
    compress_and_write_block(&[], os)
}

/// Port of `blbutil/BGZIPOutputStream.java`. Buffers up to `MAX_INPUT_BYTES` and emits
/// one BGZF block per full buffer. Call [`BgzipOutputStream::close`] to flush the final
/// (partial) block and, if requested, the empty EOF block.
pub struct BgzipOutputStream<W: Write> {
    write_empty_block: bool,
    os: W,
    input: Vec<u8>,
}

impl<W: Write> BgzipOutputStream<W> {
    /// `new BGZIPOutputStream(OutputStream os, boolean writeEmptyBlock)`.
    pub fn new(os: W, write_empty_block: bool) -> Self {
        BgzipOutputStream {
            write_empty_block,
            os,
            input: Vec::with_capacity(MAX_INPUT_BYTES),
        }
    }

    fn flush_block(&mut self) -> io::Result<()> {
        compress_and_write_block(&self.input, &mut self.os)?;
        self.input.clear();
        Ok(())
    }

    /// `close()` — flushes the final block, optionally the empty EOF block, and returns
    /// the underlying writer.
    pub fn close(mut self) -> io::Result<W> {
        if !self.input.is_empty() {
            self.flush_block()?;
        }
        if self.write_empty_block {
            self.flush_block()?; // empty buffer -> 28-byte EOF block
        }
        self.os.flush()?;
        Ok(self.os)
    }
}

impl<W: Write> Write for BgzipOutputStream<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut data = buf;
        while !data.is_empty() {
            let avail = MAX_INPUT_BYTES - self.input.len();
            let take = avail.min(data.len());
            self.input.extend_from_slice(&data[..take]);
            data = &data[take..];
            if self.input.len() == MAX_INPUT_BYTES {
                self.flush_block()?;
            }
        }
        Ok(buf.len())
    }

    /// `flush()` — compresses the current buffer (matching Java, which always emits a
    /// block on flush) and flushes the underlying writer.
    fn flush(&mut self) -> io::Result<()> {
        self.flush_block()?;
        self.os.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::MultiGzDecoder;
    use std::io::Read;

    fn deterministic_input(n: usize) -> Vec<u8> {
        (0..n).map(|i| ((i * 31 + 7) % 256) as u8).collect()
    }

    fn bgzip(input: &[u8]) -> Vec<u8> {
        let mut s = BgzipOutputStream::new(Vec::new(), true);
        s.write_all(input).unwrap();
        s.close().unwrap()
    }

    #[test]
    fn empty_block_is_canonical_bgzf_eof() {
        let mut out = Vec::new();
        write_empty_block(&mut out).unwrap();
        // The canonical 28-byte BGZF EOF marker.
        let eof: [u8; 28] = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        assert_eq!(out, eof);
    }

    #[test]
    fn roundtrip_multi_block() {
        // > 3 blocks worth, plus a partial block, to exercise block boundaries.
        let input = deterministic_input(200_000);
        let compressed = bgzip(&input);
        let mut decoded = Vec::new();
        MultiGzDecoder::new(&compressed[..])
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, input);
    }

    #[test]
    fn roundtrip_empty_and_small() {
        for n in [0usize, 1, 100, MAX_INPUT_BYTES, MAX_INPUT_BYTES + 1] {
            let input = deterministic_input(n);
            let compressed = bgzip(&input);
            let mut decoded = Vec::new();
            MultiGzDecoder::new(&compressed[..])
                .read_to_end(&mut decoded)
                .unwrap();
            assert_eq!(decoded, input, "roundtrip failed for n={n}");
        }
    }
}
