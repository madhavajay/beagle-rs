//! Faithful Rust analogues of `java.io.DataOutputStream` / `java.io.DataInputStream` — the
//! big-endian primitive read/write methods plus Java's *modified UTF-8* string format.
//!
//! Beagle's bref3 binary format is written with `DataOutputStream` and read with
//! `DataInputStream`, so byte-exact bref3 output requires reproducing these exactly:
//! - integers/shorts/longs are big-endian;
//! - `writeUTF`/`readUTF` use a 2-byte unsigned length prefix followed by *modified UTF-8*
//!   (operating on UTF-16 code units: `0x0000` and `0x0080..=0x07FF` → 2 bytes,
//!   `0x0800..=0xFFFF` → 3 bytes, otherwise 1 byte; supplementary chars become two 3-byte
//!   surrogate encodings).

use std::io::{self, Read, Write};

/// Java `DataOutputStream` write methods over any [`Write`].
pub struct DataOut<W: Write> {
    inner: W,
}

impl<W: Write> DataOut<W> {
    /// Wraps a writer.
    pub fn new(inner: W) -> Self {
        DataOut { inner }
    }

    /// Consumes the wrapper, returning the inner writer.
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Borrows the inner writer.
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// `writeByte(int)` — low 8 bits.
    pub fn write_byte(&mut self, v: i32) -> io::Result<()> {
        self.inner.write_all(&[v as u8])
    }

    /// `writeBoolean(boolean)`.
    pub fn write_boolean(&mut self, v: bool) -> io::Result<()> {
        self.inner.write_all(&[u8::from(v)])
    }

    /// `writeShort(int)` — low 16 bits, big-endian.
    pub fn write_short(&mut self, v: i32) -> io::Result<()> {
        self.inner.write_all(&(v as u16).to_be_bytes())
    }

    /// `writeInt(int)` — big-endian.
    pub fn write_int(&mut self, v: i32) -> io::Result<()> {
        self.inner.write_all(&v.to_be_bytes())
    }

    /// `writeLong(long)` — big-endian.
    pub fn write_long(&mut self, v: i64) -> io::Result<()> {
        self.inner.write_all(&v.to_be_bytes())
    }

    /// `write(byte[])`.
    pub fn write_fully(&mut self, b: &[u8]) -> io::Result<()> {
        self.inner.write_all(b)
    }

    /// `writeUTF(String)` — 2-byte unsigned length prefix + modified UTF-8 bytes.
    pub fn write_utf(&mut self, s: &str) -> io::Result<()> {
        let bytes = modified_utf8(s);
        assert!(
            bytes.len() <= 65535,
            "encoded string too long: {} bytes",
            bytes.len()
        );
        self.write_short(bytes.len() as i32)?;
        self.inner.write_all(&bytes)
    }

    /// Flushes the inner writer.
    pub fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// Encodes `s` to Java modified UTF-8 (no length prefix).
fn modified_utf8(s: &str) -> Vec<u8> {
    let mut out = Vec::with_capacity(s.len());
    for c in s.encode_utf16() {
        let c = c as u32;
        if (0x0001..=0x007F).contains(&c) {
            out.push(c as u8);
        } else if c == 0 || (0x0080..=0x07FF).contains(&c) {
            out.push((0xC0 | (c >> 6)) as u8);
            out.push((0x80 | (c & 0x3F)) as u8);
        } else {
            out.push((0xE0 | (c >> 12)) as u8);
            out.push((0x80 | ((c >> 6) & 0x3F)) as u8);
            out.push((0x80 | (c & 0x3F)) as u8);
        }
    }
    out
}

/// Java `DataInputStream` read methods over any [`Read`].
pub struct DataIn<R: Read> {
    inner: R,
}

fn eof() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "EOF")
}

impl<R: Read> DataIn<R> {
    /// Wraps a reader.
    pub fn new(inner: R) -> Self {
        DataIn { inner }
    }

    /// Borrows the inner reader.
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// `readFully(byte[] b)`.
    pub fn read_fully(&mut self, b: &mut [u8]) -> io::Result<()> {
        self.inner.read_exact(b)
    }

    /// `readByte()` — signed.
    pub fn read_byte(&mut self) -> io::Result<i8> {
        let mut b = [0u8; 1];
        self.inner.read_exact(&mut b)?;
        Ok(b[0] as i8)
    }

    /// `readUnsignedByte()`.
    pub fn read_unsigned_byte(&mut self) -> io::Result<i32> {
        let mut b = [0u8; 1];
        self.inner.read_exact(&mut b)?;
        Ok(b[0] as i32)
    }

    /// `readBoolean()`.
    pub fn read_boolean(&mut self) -> io::Result<bool> {
        Ok(self.read_unsigned_byte()? != 0)
    }

    /// `readShort()` — signed, big-endian.
    pub fn read_short(&mut self) -> io::Result<i16> {
        let mut b = [0u8; 2];
        self.inner.read_exact(&mut b)?;
        Ok(i16::from_be_bytes(b))
    }

    /// `readUnsignedShort()`.
    pub fn read_unsigned_short(&mut self) -> io::Result<i32> {
        let mut b = [0u8; 2];
        self.inner.read_exact(&mut b)?;
        Ok(u16::from_be_bytes(b) as i32)
    }

    /// `readInt()` — big-endian.
    pub fn read_int(&mut self) -> io::Result<i32> {
        let mut b = [0u8; 4];
        self.inner.read_exact(&mut b)?;
        Ok(i32::from_be_bytes(b))
    }

    /// `readLong()` — big-endian.
    pub fn read_long(&mut self) -> io::Result<i64> {
        let mut b = [0u8; 8];
        self.inner.read_exact(&mut b)?;
        Ok(i64::from_be_bytes(b))
    }

    /// `readUTF()` — 2-byte unsigned length prefix + modified UTF-8 bytes.
    pub fn read_utf(&mut self) -> io::Result<String> {
        let len = self.read_unsigned_short()? as usize;
        let mut bytes = vec![0u8; len];
        self.inner.read_exact(&mut bytes)?;
        decode_modified_utf8(&bytes)
    }
}

/// Decodes Java modified UTF-8 bytes (no length prefix) into a `String`.
fn decode_modified_utf8(bytes: &[u8]) -> io::Result<String> {
    let mut units: Vec<u16> = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        let a = bytes[i] as u32;
        if a < 0x80 {
            units.push(a as u16);
            i += 1;
        } else if (a & 0xE0) == 0xC0 {
            if i + 1 >= bytes.len() {
                return Err(eof());
            }
            let b = bytes[i + 1] as u32;
            units.push((((a & 0x1F) << 6) | (b & 0x3F)) as u16);
            i += 2;
        } else if (a & 0xF0) == 0xE0 {
            if i + 2 >= bytes.len() {
                return Err(eof());
            }
            let b = bytes[i + 1] as u32;
            let c = bytes[i + 2] as u32;
            units.push((((a & 0x0F) << 12) | ((b & 0x3F) << 6) | (c & 0x3F)) as u16);
            i += 3;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "malformed modified UTF-8",
            ));
        }
    }
    String::from_utf16(&units)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid UTF-16"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn primitive_round_trip() {
        let mut out = DataOut::new(Vec::new());
        out.write_int(0x0102_0304).unwrap();
        out.write_int(-1).unwrap();
        out.write_short(0x0506).unwrap();
        out.write_byte(0xFF).unwrap();
        out.write_long(0x0102_0304_0506_0708).unwrap();
        out.write_boolean(true).unwrap();
        let bytes = out.into_inner();
        // big-endian layout
        assert_eq!(&bytes[0..4], &[0x01, 0x02, 0x03, 0x04]);
        assert_eq!(&bytes[4..8], &[0xFF, 0xFF, 0xFF, 0xFF]);
        assert_eq!(&bytes[8..10], &[0x05, 0x06]);
        assert_eq!(bytes[10], 0xFF);

        let mut din = DataIn::new(Cursor::new(bytes));
        assert_eq!(din.read_int().unwrap(), 0x0102_0304);
        assert_eq!(din.read_int().unwrap(), -1);
        assert_eq!(din.read_unsigned_short().unwrap(), 0x0506);
        assert_eq!(din.read_unsigned_byte().unwrap(), 0xFF);
        assert_eq!(din.read_long().unwrap(), 0x0102_0304_0506_0708);
        assert!(din.read_boolean().unwrap());
    }

    #[test]
    fn write_utf_known_bytes() {
        // "AB" -> length 2, bytes 'A','B'
        let mut out = DataOut::new(Vec::new());
        out.write_utf("AB").unwrap();
        assert_eq!(out.into_inner(), vec![0x00, 0x02, b'A', b'B']);

        // U+00E9 (é) is in 0x0080..=0x07FF -> 2 bytes: 0xC3 0xA9; length prefix 0x0002
        let mut out = DataOut::new(Vec::new());
        out.write_utf("\u{00E9}").unwrap();
        assert_eq!(out.into_inner(), vec![0x00, 0x02, 0xC3, 0xA9]);

        // U+0000 (NUL) encodes as 2 bytes 0xC0 0x80 in modified UTF-8
        let mut out = DataOut::new(Vec::new());
        out.write_utf("\u{0000}").unwrap();
        assert_eq!(out.into_inner(), vec![0x00, 0x02, 0xC0, 0x80]);
    }

    #[test]
    fn utf_round_trip_incl_3byte_and_surrogates() {
        let samples = ["chr1", "rs12345", "é-ñ", "\u{0000}embedded", "𝔘nicode"];
        for s in samples {
            let mut out = DataOut::new(Vec::new());
            out.write_utf(s).unwrap();
            let mut din = DataIn::new(Cursor::new(out.into_inner()));
            assert_eq!(din.read_utf().unwrap(), s);
        }
    }

    #[test]
    fn read_fully_and_byte_signedness() {
        let mut din = DataIn::new(Cursor::new(vec![0x80u8, 0x7F, 0xFE]));
        assert_eq!(din.read_byte().unwrap(), -128i8);
        let mut buf = [0u8; 2];
        din.read_fully(&mut buf).unwrap();
        assert_eq!(buf, [0x7F, 0xFE]);
    }
}
