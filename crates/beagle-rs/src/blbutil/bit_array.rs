//! Port of `blbutil/BitArray.java` — a fixed-length mutable bit sequence backed by
//! `long[]` words. All shifts mirror Java `long` semantics (shift count masked to 6
//! bits) so the word layout is byte-exact, which matters for bref3 serialization.

use super::{jshl_l, jshr_l, jushr_l};

const LOG2_BITS_PER_WORD: i32 = 6;
const WORD_MASK: i64 = -1; // 0xffffffffffffffff
const BIT_INDEX_MASK: i32 = (1 << LOG2_BITS_PER_WORD) - 1;

/// Port of `blbutil/BitArray.java`.
#[derive(Clone)]
pub struct BitArray {
    words: Vec<i64>,
    size: i32,
}

#[inline]
fn n_words(size: i32) -> usize {
    // (size + Long.SIZE - 1) / Long.SIZE
    ((size + 64 - 1) / 64) as usize
}

impl BitArray {
    /// `new BitArray(int size)` — all bits 0.
    pub fn new(size: i32) -> Self {
        assert!(size >= 0, "{}", size);
        BitArray {
            words: vec![0i64; n_words(size)],
            size,
        }
    }

    /// `new BitArray(long[] values, int size)`.
    pub fn from_words(values: &[i64], size: i32) -> Self {
        assert!(size >= 0, "{}", size);
        let nw = n_words(size);
        assert!(values.len() == nw, "{}", values.len());
        BitArray {
            words: values.to_vec(),
            size,
        }
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.size
    }

    /// `get(int index)`.
    pub fn get(&self, index: i32) -> bool {
        assert!(index < self.size, "{}", index);
        let wi = (index >> LOG2_BITS_PER_WORD) as usize;
        (self.words[wi] & jshl_l(1, index)) != 0
    }

    /// `getAsInt(int index)` — note: for `index % 64 == 63` a set bit yields `-1`
    /// (the arithmetic `>>` of the sign bit), faithfully preserved.
    pub fn get_as_int(&self, index: i32) -> i32 {
        assert!(index < self.size, "{}", index);
        let wi = (index >> LOG2_BITS_PER_WORD) as usize;
        jshr_l(self.words[wi] & jshl_l(1, index), index) as i32
    }

    /// `set(int index)`.
    pub fn set(&mut self, index: i32) {
        assert!(index < self.size, "{}", index);
        let wi = (index >> LOG2_BITS_PER_WORD) as usize;
        self.words[wi] |= jshl_l(1, index);
    }

    /// `clear(int index)`.
    pub fn clear(&mut self, index: i32) {
        assert!(index < self.size, "{}", index);
        let wi = (index >> LOG2_BITS_PER_WORD) as usize;
        self.words[wi] &= !jshl_l(1, index);
    }

    /// `clear()` — clears all bits.
    pub fn clear_all(&mut self) {
        self.words.iter_mut().for_each(|w| *w = 0);
    }

    /// `restrict(int from, int to)` — a new `BitArray` of the bits `[from, to)`.
    pub fn restrict(&self, from: i32, to: i32) -> BitArray {
        assert!(!(from < 0 || from > to || to > self.size), "{}", from);
        if from == to {
            return BitArray::new(0);
        }
        let mut result = BitArray::new(to - from);
        let n_result_words_m1 = result.words.len() - 1;
        let is_word_aligned = (from & BIT_INDEX_MASK) == 0;
        let mut src_word = (from >> LOG2_BITS_PER_WORD) as usize;
        for w in 0..n_result_words_m1 {
            result.words[w] = if is_word_aligned {
                self.words[src_word]
            } else {
                jushr_l(self.words[src_word], from) | jshl_l(self.words[src_word + 1], -from)
            };
            src_word += 1;
        }
        let end_word_mask = jushr_l(WORD_MASK, -to);
        result.words[n_result_words_m1] = if ((to - 1) & BIT_INDEX_MASK) < (from & BIT_INDEX_MASK) {
            jushr_l(self.words[src_word], from)
                | jshl_l(self.words[src_word + 1] & end_word_mask, -from)
        } else {
            jushr_l(self.words[src_word] & end_word_mask, from)
        };
        result
    }

    /// `copyFrom(BitArray src, int from, int to)`.
    pub fn copy_from(&mut self, src: &BitArray, from: i32, to: i32) {
        assert!(
            !(from < 0 || from > to || to > self.size || to > src.size()),
            "{}",
            from
        );
        if from == to {
            return;
        }
        let start_word = (from >> LOG2_BITS_PER_WORD) as usize;
        let end_word = ((to - 1) >> LOG2_BITS_PER_WORD) as usize;
        let start_word_mask = jshl_l(WORD_MASK, from);
        let end_word_mask = jushr_l(WORD_MASK, -to);
        if start_word == end_word {
            let mask = start_word_mask & end_word_mask;
            self.words[start_word] ^= (self.words[start_word] ^ src.words[start_word]) & mask;
        } else {
            self.words[start_word] ^=
                (self.words[start_word] ^ src.words[start_word]) & start_word_mask;
            for j in (start_word + 1)..end_word {
                self.words[j] = src.words[j];
            }
            self.words[end_word] ^= (self.words[end_word] ^ src.words[end_word]) & end_word_mask;
        }
    }

    /// `hash(int from, int to)`.
    pub fn hash(&self, from: i32, to: i32) -> i32 {
        assert!(!(from < 0 || from > to || to > self.size), "{}", from);
        if from == to {
            return 0;
        }
        let start_word = (from >> LOG2_BITS_PER_WORD) as usize;
        let end_word = ((to - 1) >> LOG2_BITS_PER_WORD) as usize;
        let start_word_mask = jshl_l(WORD_MASK, from);
        let end_word_mask = jushr_l(WORD_MASK, -to);
        if start_word == end_word {
            let mask = start_word_mask & end_word_mask;
            Self::long_hash_code(self.words[start_word] & mask)
        } else {
            let mut long_hash = self.words[start_word] & start_word_mask;
            for j in (start_word + 1)..end_word {
                long_hash ^= self.words[j];
            }
            long_hash ^= self.words[end_word] & end_word_mask;
            Self::long_hash_code(long_hash)
        }
    }

    /// `longHashCode(long value)`.
    pub fn long_hash_code(value: i64) -> i32 {
        (value ^ jushr_l(value, 32)) as i32
    }

    /// `swapBits(BitArray a, BitArray b, int from, int to)`.
    pub fn swap_bits(a: &mut BitArray, b: &mut BitArray, from: i32, to: i32) {
        assert!(a.size() == b.size(), "inconsistent data");
        assert!(!(from < 0 || from > to || to > a.size()), "{}", from);
        if from == to {
            return;
        }
        let start_word = (from >> LOG2_BITS_PER_WORD) as usize;
        let end_word = ((to - 1) >> LOG2_BITS_PER_WORD) as usize;
        let start_word_mask = jshl_l(WORD_MASK, from);
        let end_word_mask = jushr_l(WORD_MASK, -to);
        if start_word == end_word {
            let mask = start_word_mask & end_word_mask;
            a.words[start_word] ^= b.words[start_word] & mask;
            b.words[start_word] ^= a.words[start_word] & mask;
            a.words[start_word] ^= b.words[start_word] & mask;
        } else {
            a.words[start_word] ^= b.words[start_word] & start_word_mask;
            b.words[start_word] ^= a.words[start_word] & start_word_mask;
            a.words[start_word] ^= b.words[start_word] & start_word_mask;
            for j in (start_word + 1)..end_word {
                // full-word xor-swap == swapping the words
                std::mem::swap(&mut a.words[j], &mut b.words[j]);
            }
            a.words[end_word] ^= b.words[end_word] & end_word_mask;
            b.words[end_word] ^= a.words[end_word] & end_word_mask;
            a.words[end_word] ^= b.words[end_word] & end_word_mask;
        }
    }

    /// `equal(BitArray other, int from, int to)`.
    pub fn equal(&self, other: &BitArray, from: i32, to: i32) -> bool {
        assert!(
            !(from < 0 || from > to || to > self.size || to > other.size()),
            "{}",
            from
        );
        if from == to {
            return true;
        }
        let start_word = (from >> LOG2_BITS_PER_WORD) as usize;
        let end_word = ((to - 1) >> LOG2_BITS_PER_WORD) as usize;
        let start_word_mask = jshl_l(WORD_MASK, from);
        let end_word_mask = jushr_l(WORD_MASK, -to);
        if start_word == end_word {
            let mask = start_word_mask & end_word_mask;
            ((self.words[start_word] ^ other.words[start_word]) & mask) == 0
        } else {
            let mut are_equal =
                ((self.words[start_word] ^ other.words[start_word]) & start_word_mask) == 0;
            for j in (start_word + 1)..end_word {
                are_equal &= self.words[j] == other.words[j];
            }
            are_equal &= ((self.words[end_word] ^ other.words[end_word]) & end_word_mask) == 0;
            are_equal
        }
    }

    /// `toLongArray()`.
    pub fn to_long_array(&self) -> Vec<i64> {
        self.words.clone()
    }

    /// `equals(BitArray a, BitArray b)`.
    pub fn equals(a: &BitArray, b: &BitArray) -> bool {
        a.size() == b.size() && a.words == b.words
    }
}

impl std::fmt::Display for BitArray {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::with_capacity(self.size.max(0) as usize);
        for j in 0..self.size {
            s.push(if self.get(j) { '1' } else { '0' });
        }
        f.write_str(&s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_get_clear() {
        let mut b = BitArray::new(200);
        assert_eq!(b.size(), 200);
        b.set(0);
        b.set(63);
        b.set(64);
        b.set(199);
        assert!(b.get(0) && b.get(63) && b.get(64) && b.get(199));
        assert!(!b.get(1) && !b.get(62) && !b.get(65));
        assert_eq!(b.get_as_int(0), 1);
        assert_eq!(b.get_as_int(64), 1);
        // ints quirk: bit 63 (index % 64 == 63) yields -1 via arithmetic shift
        assert_eq!(b.get_as_int(63), -1);
        b.clear(63);
        assert!(!b.get(63));
        b.clear_all();
        assert!(!b.get(0) && !b.get(199));
    }

    #[test]
    fn restrict_matches_manual_bits() {
        // Build a 130-bit array with a known pattern, restrict, and compare bit-by-bit.
        let mut b = BitArray::new(130);
        for j in (0..130).filter(|j| j % 3 == 0) {
            b.set(j);
        }
        let r = b.restrict(5, 128);
        assert_eq!(r.size(), 123);
        for j in 0..r.size() {
            assert_eq!(r.get(j), b.get(j + 5), "bit {j}");
        }
    }

    #[test]
    fn copy_from_and_equal() {
        let mut src = BitArray::new(128);
        for j in (0..128).filter(|j| j % 5 == 0) {
            src.set(j);
        }
        let mut dst = BitArray::new(128);
        dst.copy_from(&src, 10, 100);
        assert!(dst.equal(&src, 10, 100));
        assert!(!dst.get(0)); // outside the copied range stays 0
    }

    #[test]
    fn swap_bits_roundtrip() {
        let mut a = BitArray::new(128);
        let mut b = BitArray::new(128);
        for j in 0..128 {
            if j % 2 == 0 {
                a.set(j);
            } else {
                b.set(j);
            }
        }
        let a0 = a.to_long_array();
        let b0 = b.to_long_array();
        BitArray::swap_bits(&mut a, &mut b, 0, 128);
        assert_eq!(a.to_long_array(), b0);
        assert_eq!(b.to_long_array(), a0);
    }

    #[test]
    fn hash_and_equals() {
        let mut a = BitArray::new(70);
        let mut b = BitArray::new(70);
        a.set(3);
        a.set(65);
        b.set(3);
        b.set(65);
        assert!(BitArray::equals(&a, &b));
        assert_eq!(a.hash(0, 70), b.hash(0, 70));
        b.set(10);
        assert!(!BitArray::equals(&a, &b));
    }

    #[test]
    fn display_bits() {
        let mut b = BitArray::new(5);
        b.set(0);
        b.set(2);
        b.set(4);
        assert_eq!(b.to_string(), "10101");
    }
}
