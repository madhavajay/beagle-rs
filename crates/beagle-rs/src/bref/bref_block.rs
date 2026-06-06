//! Port of `bref/BrefBlock.java` — the start chromosome coordinates and file offset of a
//! bref data block.

use crate::beagleutil::ChromIds;
use crate::blbutil::consts;
use std::fmt;

/// Port of `bref/BrefBlock.java`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BrefBlock {
    chrom_index: i32,
    pos: i32,
    offset: i64,
}

impl BrefBlock {
    /// `new BrefBlock(int chromIndex, int pos, long offset)`.
    pub fn new(chrom_index: i32, pos: i32, offset: i64) -> Self {
        BrefBlock {
            chrom_index,
            pos,
            offset,
        }
    }

    /// `chromIndex()`.
    pub fn chrom_index(&self) -> i32 {
        self.chrom_index
    }

    /// `pos()`.
    pub fn pos(&self) -> i32 {
        self.pos
    }

    /// `offset()`.
    pub fn offset(&self) -> i64 {
        self.offset
    }
}

impl fmt::Display for BrefBlock {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}{}{}{}{}]",
            ChromIds::instance().id(self.chrom_index),
            consts::TAB,
            self.pos,
            consts::TAB,
            self.offset
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accessors_and_display() {
        let idx = ChromIds::instance().get_index("chrBREFBLK");
        let b = BrefBlock::new(idx, 12345, 9_000_000_000);
        assert_eq!(b.chrom_index(), idx);
        assert_eq!(b.pos(), 12345);
        assert_eq!(b.offset(), 9_000_000_000);
        assert_eq!(b.to_string(), "[chrBREFBLK\t12345\t9000000000]");
    }
}
