//! Port of `blbutil/BlockLineReader.java` — reads a file as successive blocks of up to
//! `blockSize` lines, preserving order.
//!
//! Java prefetches blocks on a background thread into a bounded queue; that only affects
//! timing, never the block *sequence*, so this port reads blocks sequentially. An empty
//! block (`Vec::is_empty()`) is the end-of-input sentinel (Java's unique length-0 array).

use std::path::Path;

use super::FileIt;

/// Port of `blbutil/BlockLineReader.java`.
pub struct BlockLineReader<I: FileIt<Item = String>> {
    it: I,
    block_size: i32,
}

impl<I: FileIt<Item = String>> BlockLineReader<I> {
    /// `BlockLineReader.create(FileIt<String> it, int blockSize, int nBlocks)`. The
    /// `n_blocks` prefetch-depth hint does not affect output.
    pub fn create(it: I, block_size: i32, n_blocks: i32) -> Self {
        assert!(block_size >= 1, "{block_size}");
        assert!(n_blocks >= 1, "{n_blocks}");
        BlockLineReader { it, block_size }
    }

    /// `next()` — the next block of up to `blockSize` lines, or an empty `Vec` (the
    /// sentinel) once the input is exhausted.
    pub fn next_block(&mut self) -> Vec<String> {
        let mut buffer = Vec::with_capacity(self.block_size as usize);
        while (buffer.len() as i32) < self.block_size {
            match self.it.next() {
                Some(line) => buffer.push(line),
                None => break,
            }
        }
        buffer
    }

    /// `file()`.
    pub fn file(&self) -> Option<&Path> {
        self.it.file()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blbutil::InputIt;
    use std::io::{BufReader, Cursor};

    fn reader(text: &str, block_size: i32) -> BlockLineReader<InputIt> {
        let it = InputIt::from_reader(Box::new(BufReader::new(Cursor::new(text.to_owned()))), None);
        BlockLineReader::create(it, block_size, 1)
    }

    #[test]
    fn yields_blocks_then_sentinel() {
        let mut r = reader("a\nb\nc\nd\ne\n", 2);
        assert_eq!(r.next_block(), vec!["a", "b"]);
        assert_eq!(r.next_block(), vec!["c", "d"]);
        assert_eq!(r.next_block(), vec!["e"]); // partial final block
        assert!(r.next_block().is_empty()); // sentinel
        assert!(r.next_block().is_empty()); // sentinel repeats
    }

    #[test]
    fn block_size_larger_than_input() {
        let mut r = reader("x\ny\n", 10);
        assert_eq!(r.next_block(), vec!["x", "y"]);
        assert!(r.next_block().is_empty());
    }
}
