//! Port of the `IntArray` implementors in the `ints` package:
//! `PackedIntArray`, `UnsignedByteArray`, `CharArray`, `WrappedIntArray`, `IndexArray`.

use super::{
    as_string, int_slice_to_string, java_binary_search, jshl, jshr, jushr, IntArray, IntList,
};

/// `java.util.Arrays.copyOfRange(byte[], from, to)` semantics (pad with 0 if `to`
/// exceeds the source length).
fn byte_copy_of_range(src: &[u8], from: i32, to: i32) -> Vec<u8> {
    let new_len = to - from; // IllegalArgumentException analog if to < from
    assert!(from >= 0 && from as usize <= src.len(), "{}", from);
    let mut out = vec![0u8; new_len as usize];
    let n = ((src.len() as i32 - from).min(new_len)).max(0) as usize;
    out[0..n].copy_from_slice(&src[from as usize..from as usize + n]);
    out
}

// =================================================================================
// PackedIntArray — ints/PackedIntArray.java
// =================================================================================

/// Port of `ints/PackedIntArray.java`: an immutable array of non-negative ints
/// packed into 1/2/4/… bits per value inside a backing `i32[]`.
pub struct PackedIntArray {
    pack_index: i32,  // each value is stored in (1 << pack_index) bits
    index_shift: i32, // right-shift mapping array index → packed-int index
    values_per_int_m1: i32,
    size: i32,
    ia: Vec<i32>,
}

/// `numberOfTrailingZeros(Integer.SIZE)` = `numberOfTrailingZeros(32)` = 5.
const MAX_PACK_INDEX: i32 = 5;

impl PackedIntArray {
    fn round_up_to_power_of_two(mut x: i32) -> i32 {
        x = x.wrapping_sub(1);
        x |= jshr(x, 1);
        x |= jshr(x, 2);
        x |= jshr(x, 4);
        x |= jshr(x, 8);
        x |= jshr(x, 16);
        x.wrapping_add(1)
    }

    /// `packIndex(valueSize)` — log2 of the (power-of-two) bit width per value.
    fn pack_index(value_size: i32) -> i32 {
        assert!(value_size >= 1, "{}", value_size);
        if value_size == 1 {
            0
        } else {
            let next = Self::round_up_to_power_of_two(value_size);
            let n_mask_bits = next.trailing_zeros() as i32;
            let next2 = Self::round_up_to_power_of_two(n_mask_bits);
            next2.trailing_zeros() as i32
        }
    }

    fn with_fields(ia: Vec<i32>, size: i32, pack_index: i32) -> Self {
        PackedIntArray {
            pack_index,
            index_shift: MAX_PACK_INDEX - pack_index,
            values_per_int_m1: jshr(32, pack_index) - 1,
            size,
            ia,
        }
    }

    /// `new PackedIntArray(int[] ia, int valueSize)`.
    pub fn from_slice(values: &[i32], value_size: i32) -> Self {
        Self::from_iter_checked(values.iter().copied(), values.len() as i32, value_size)
    }

    /// `new PackedIntArray(IntList il, int valueSize)`.
    pub fn from_list(il: &IntList, value_size: i32) -> Self {
        Self::from_iter_checked((0..il.size()).map(|j| il.get(j)), il.size(), value_size)
    }

    fn from_iter_checked<I: Iterator<Item = i32>>(values: I, size: i32, value_size: i32) -> Self {
        assert!(value_size >= 1, "{}", value_size);
        let pack_index = Self::pack_index(value_size);
        let index_shift = MAX_PACK_INDEX - pack_index;
        let values_per_int_m1 = jshr(32, pack_index) - 1;
        let words = ((size + values_per_int_m1) / (values_per_int_m1 + 1)) as usize;
        let mut ia = vec![0i32; words];
        for (j, value) in values.enumerate() {
            let j = j as i32;
            assert!(value >= 0 && value < value_size, "{}", value);
            let word = jshr(j, index_shift) as usize;
            ia[word] |= jshl(value, jshl(j & values_per_int_m1, pack_index));
        }
        Self::with_fields(ia, size, pack_index)
    }

    /// `fromSignedByteArray(byte[] ba, int valueSize)`. Note: despite the name, the
    /// "signed" path masks each byte with `0xff` (ints quirk — preserved).
    pub fn from_signed_byte_array(ba: &[u8], value_size: i32) -> Self {
        Self::from_byte_array(ba, 0, ba.len() as i32, value_size, false)
    }

    /// `fromSignedByteArray(byte[] ba, int from, int to, int valueSize)`.
    pub fn from_signed_byte_array_range(ba: &[u8], from: i32, to: i32, value_size: i32) -> Self {
        Self::from_byte_array(ba, from, to, value_size, false)
    }

    /// `fromUnsignedByteArray(byte[] ba, int valueSize)`. Note: the "unsigned" path
    /// masks each byte with `Byte.MAX_VALUE` (`0x7f`) — preserved quirk.
    pub fn from_unsigned_byte_array(ba: &[u8], value_size: i32) -> Self {
        Self::from_byte_array(ba, 0, ba.len() as i32, value_size, true)
    }

    /// `fromUnsignedByteArray(byte[] ba, int from, int to, int valueSize)`.
    pub fn from_unsigned_byte_array_range(ba: &[u8], from: i32, to: i32, value_size: i32) -> Self {
        Self::from_byte_array(ba, from, to, value_size, true)
    }

    fn from_byte_array(ba: &[u8], from: i32, to: i32, value_size: i32, use_unsigned: bool) -> Self {
        assert!(value_size >= 1, "{}", value_size);
        let mask: i32 = if use_unsigned { i8::MAX as i32 } else { 0xff };
        let pack_index = Self::pack_index(value_size);
        let values_per_int_m1 = jshr(32, pack_index) - 1;
        let index_shift = MAX_PACK_INDEX - pack_index;

        let size = to - from;
        let words = ((size + values_per_int_m1) / (values_per_int_m1 + 1)) as usize;
        let mut ia = vec![0i32; words];
        for j in from..to {
            let offset = j - from;
            // Java byte is signed: replicate sign-extension before masking.
            let value = ((ba[j as usize] as i8) as i32) & mask;
            assert!(value >= 0 && value < value_size, "{}", value);
            // ints-6: word index uses absolute `j`, bit offset uses relative `offset`.
            let word = jshr(j, index_shift) as usize;
            ia[word] |= jshl(value, jshl(offset & values_per_int_m1, pack_index));
        }
        Self::with_fields(ia, size, pack_index)
    }

    /// `fromUnsignedTwoByteArray(byte[] ba, int valueSize)` — big-endian 2-byte values.
    pub fn from_unsigned_two_byte_array(ba: &[u8], value_size: i32) -> Self {
        assert!(value_size >= 1, "{}", value_size);
        assert!((ba.len() & 1) == 0, "{}", ba.len());
        let pack_index = Self::pack_index(value_size);
        let index_shift = MAX_PACK_INDEX - pack_index;
        let values_per_int_m1 = jshr(32, pack_index) - 1;

        let size = (ba.len() / 2) as i32;
        let words = ((size + values_per_int_m1) / (values_per_int_m1 + 1)) as usize;
        let mut ia = vec![0i32; words];
        for j in 0..size {
            let hi = (ba[(2 * j) as usize] as i32) & 0xff;
            let lo = (ba[(2 * j + 1) as usize] as i32) & 0xff;
            let value = 0xffff & ((hi << 8) | lo);
            assert!(value >= 0 && value < value_size, "{}", value);
            let word = jshr(j, index_shift) as usize;
            ia[word] |= jshl(value, jshl(j & values_per_int_m1, pack_index));
        }
        Self::with_fields(ia, size, pack_index)
    }

    /// The backing packed `int[]` words (for serialization and parity checks).
    pub fn words(&self) -> &[i32] {
        &self.ia
    }
}

impl IntArray for PackedIntArray {
    fn size(&self) -> i32 {
        self.size
    }

    fn get(&self, index: i32) -> i32 {
        // ints-1: Java checks `index > size` (not `>=`) and ignores negatives.
        assert!(index <= self.size, "{}", index);
        let bits_per_value = jshl(1, self.pack_index);
        let value_mask = jshl(1, bits_per_value) - 1;
        let word = self.ia[jshr(index, self.index_shift) as usize];
        value_mask & jushr(word, jshl(index & self.values_per_int_m1, self.pack_index))
    }
}

// =================================================================================
// UnsignedByteArray — ints/UnsignedByteArray.java
// =================================================================================

/// Port of `ints/UnsignedByteArray.java`: ints in `0..=255` stored as bytes.
pub struct UnsignedByteArray {
    ba: Vec<u8>,
}

impl UnsignedByteArray {
    /// `new UnsignedByteArray(byte[] ba)` — clones, interpreting bytes as unsigned.
    pub fn from_bytes(ba: &[u8]) -> Self {
        UnsignedByteArray { ba: ba.to_vec() }
    }

    /// `new UnsignedByteArray(byte[] ba, int from, int to)`.
    pub fn from_bytes_range(ba: &[u8], from: i32, to: i32) -> Self {
        UnsignedByteArray {
            ba: byte_copy_of_range(ba, from, to),
        }
    }

    /// `new UnsignedByteArray(int[] ia)` — values must be in `0..=255`.
    pub fn from_ints(ia: &[i32]) -> Self {
        Self::from_ints_range(ia, 0, ia.len() as i32)
    }

    /// `new UnsignedByteArray(IntList il)`.
    pub fn from_list(il: &IntList) -> Self {
        Self::from_list_range(il, 0, il.size())
    }

    /// `new UnsignedByteArray(int[] ia, int valueSize)`.
    pub fn from_ints_value_size(ia: &[i32], value_size: i32) -> Self {
        assert!(value_size >= 1 && value_size <= 256, "{}", value_size);
        let mut ba = vec![0u8; ia.len()];
        for (j, &v) in ia.iter().enumerate() {
            assert!(v >= 0 && v < value_size, "{}", v);
            ba[j] = v as u8;
        }
        UnsignedByteArray { ba }
    }

    /// `new UnsignedByteArray(IntList il, int valueSize)`.
    pub fn from_list_value_size(il: &IntList, value_size: i32) -> Self {
        assert!(value_size >= 1 && value_size <= 256, "{}", value_size);
        let mut ba = vec![0u8; il.size() as usize];
        for j in 0..il.size() {
            let v = il.get(j);
            assert!(v >= 0 && v < value_size, "{}", v);
            ba[j as usize] = v as u8;
        }
        UnsignedByteArray { ba }
    }

    /// `new UnsignedByteArray(int[] ia, int from, int to)` — values in `0..=255`.
    pub fn from_ints_range(ia: &[i32], from: i32, to: i32) -> Self {
        let mut ba = vec![0u8; (to - from) as usize];
        for j in from..to {
            let v = ia[j as usize];
            assert!(v >= 0 && v <= 255, "{}", v);
            ba[(j - from) as usize] = v as u8;
        }
        UnsignedByteArray { ba }
    }

    /// `new UnsignedByteArray(IntList il, int from, int to)`.
    pub fn from_list_range(il: &IntList, from: i32, to: i32) -> Self {
        let mut ba = vec![0u8; (to - from) as usize];
        for j in from..to {
            let v = il.get(j);
            assert!(v >= 0 && v <= 255, "{}", v);
            ba[(j - from) as usize] = v as u8;
        }
        UnsignedByteArray { ba }
    }

    /// `new UnsignedByteArray(ByteArrayOutputStream baos)`.
    pub fn from_baos(bytes: Vec<u8>) -> Self {
        UnsignedByteArray { ba: bytes }
    }

    /// The backing bytes (for `write(OutputStream)`).
    pub fn as_bytes(&self) -> &[u8] {
        &self.ba
    }
}

impl IntArray for UnsignedByteArray {
    fn size(&self) -> i32 {
        self.ba.len() as i32
    }

    fn get(&self, index: i32) -> i32 {
        self.ba[index as usize] as i32
    }
}

impl std::fmt::Display for UnsignedByteArray {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&as_string(self))
    }
}

// =================================================================================
// CharArray — ints/CharArray.java
// =================================================================================

/// Port of `ints/CharArray.java`: ints in `0..=65535` stored as `char` (u16).
pub struct CharArray {
    ca: Vec<u16>,
}

/// `Character.MAX_VALUE`.
const CHAR_MAX_VALUE: i32 = 65535;

impl CharArray {
    /// `new CharArray(byte[] ba)` — big-endian char pairs. ints-2: the `% 1` length
    /// check is dead code, so odd-length input panics on the trailing pair read.
    pub fn from_bytes(ba: &[u8]) -> Self {
        let mut ca = vec![0u16; ba.len() / 2];
        let mut j = 0usize;
        while j < ba.len() {
            let b1 = (ba[j] as i32) & 0xff;
            let b2 = (ba[j + 1] as i32) & 0xff;
            ca[j / 2] = ((b1 << 8) + b2) as u16;
            j += 2;
        }
        CharArray { ca }
    }

    /// `new CharArray(char[] ca)`.
    pub fn from_chars(ca: &[u16]) -> Self {
        CharArray { ca: ca.to_vec() }
    }

    /// `new CharArray(int[] ia)` — values in `0..=65535`.
    pub fn from_ints(ia: &[i32]) -> Self {
        Self::from_ints_range(ia, 0, ia.len() as i32)
    }

    /// `new CharArray(IntList il)`.
    pub fn from_list(il: &IntList) -> Self {
        Self::from_list_range(il, 0, il.size())
    }

    /// `new CharArray(int[] ia, int valueSize)`.
    pub fn from_ints_value_size(ia: &[i32], value_size: i32) -> Self {
        assert!(
            value_size >= 1 && value_size <= CHAR_MAX_VALUE + 1,
            "{}",
            value_size
        );
        let mut ca = vec![0u16; ia.len()];
        for (j, &v) in ia.iter().enumerate() {
            assert!(v >= 0 && v < value_size, "{}", v);
            ca[j] = v as u16;
        }
        CharArray { ca }
    }

    /// `new CharArray(IntList il, int valueSize)`.
    pub fn from_list_value_size(il: &IntList, value_size: i32) -> Self {
        assert!(
            value_size >= 1 && value_size <= CHAR_MAX_VALUE + 1,
            "{}",
            value_size
        );
        let mut ca = vec![0u16; il.size() as usize];
        for j in 0..il.size() {
            let v = il.get(j);
            assert!(v >= 0 && v < value_size, "{}", v);
            ca[j as usize] = v as u16;
        }
        CharArray { ca }
    }

    /// `new CharArray(int[] ia, int to, int from)` — ints-3: positionally `(start, end)`
    /// with confusingly swapped parameter names. Values in `0..=65535`.
    pub fn from_ints_range(ia: &[i32], start: i32, end: i32) -> Self {
        let mut ca = vec![0u16; (end - start) as usize];
        for j in start..end {
            let v = ia[j as usize];
            assert!(v >= 0 && v <= CHAR_MAX_VALUE, "{}", v);
            ca[(j - start) as usize] = v as u16;
        }
        CharArray { ca }
    }

    /// `new CharArray(IntList il, int from, int to)`.
    pub fn from_list_range(il: &IntList, from: i32, to: i32) -> Self {
        let mut ca = vec![0u16; (to - from) as usize];
        for j in from..to {
            let v = il.get(j);
            assert!(v >= 0 && v <= CHAR_MAX_VALUE, "{}", v);
            ca[(j - from) as usize] = v as u16;
        }
        CharArray { ca }
    }
}

impl IntArray for CharArray {
    fn size(&self) -> i32 {
        self.ca.len() as i32
    }

    fn get(&self, index: i32) -> i32 {
        self.ca[index as usize] as i32
    }
}

// =================================================================================
// WrappedIntArray — ints/WrappedIntArray.java
// =================================================================================

/// Port of `ints/WrappedIntArray.java`: an immutable `int[]` (no packing).
pub struct WrappedIntArray {
    ia: Vec<i32>,
}

impl WrappedIntArray {
    /// `new WrappedIntArray(int[] ia)`.
    pub fn from_slice(ia: &[i32]) -> Self {
        WrappedIntArray { ia: ia.to_vec() }
    }

    /// `new WrappedIntArray(IntStream values)`.
    pub fn from_int_stream<I: IntoIterator<Item = i32>>(values: I) -> Self {
        WrappedIntArray {
            ia: values.into_iter().collect(),
        }
    }

    /// `new WrappedIntArray(int[] ia, int valueSize)`.
    pub fn from_slice_value_size(ia: &[i32], value_size: i32) -> Self {
        let mut out = vec![0i32; ia.len()];
        for (j, &v) in ia.iter().enumerate() {
            assert!(v >= 0 && v < value_size, "{}", v);
            out[j] = v;
        }
        WrappedIntArray { ia: out }
    }

    /// `new WrappedIntArray(IntList il)`.
    pub fn from_list(il: &IntList) -> Self {
        WrappedIntArray { ia: il.to_array() }
    }

    /// `new WrappedIntArray(IntList il, int valueSize)`.
    pub fn from_list_value_size(il: &IntList, value_size: i32) -> Self {
        let mut out = vec![0i32; il.size() as usize];
        for j in 0..il.size() {
            let v = il.get(j);
            assert!(v >= 0 && v < value_size, "{}", v);
            out[j as usize] = v;
        }
        WrappedIntArray { ia: out }
    }

    /// `toArray()`.
    pub fn to_array(&self) -> Vec<i32> {
        self.ia.clone()
    }

    /// `binarySearch(int key)` — Java semantics over the whole array.
    pub fn binary_search(&self, key: i32) -> i32 {
        java_binary_search(&self.ia, 0, self.ia.len() as i32, key)
    }

    /// `binarySearch(int fromIndex, int toIndex, int key)`.
    pub fn binary_search_range(&self, from_index: i32, to_index: i32, key: i32) -> i32 {
        assert!(
            from_index <= to_index,
            "fromIndex({}) > toIndex({})",
            from_index,
            to_index
        );
        java_binary_search(&self.ia, from_index, to_index, key)
    }
}

impl IntArray for WrappedIntArray {
    fn size(&self) -> i32 {
        self.ia.len() as i32
    }

    fn get(&self, index: i32) -> i32 {
        self.ia[index as usize]
    }
}

impl std::fmt::Display for WrappedIntArray {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&int_slice_to_string(&self.ia))
    }
}

// =================================================================================
// IndexArray — ints/IndexArray.java
// =================================================================================

/// Port of `ints/IndexArray.java`: an [`IntArray`] paired with a value-size bound.
pub struct IndexArray {
    int_array: Box<dyn IntArray>,
    value_size: i32,
}

impl IndexArray {
    /// `new IndexArray(int[] intArray, int valueSize)`.
    pub fn from_ints(int_array: &[i32], value_size: i32) -> Self {
        IndexArray {
            int_array: super::packed_create(int_array, value_size),
            value_size,
        }
    }

    /// `new IndexArray(IntArray intArray, int valueSize)`.
    pub fn from_int_array(int_array: Box<dyn IntArray>, value_size: i32) -> Self {
        IndexArray {
            int_array,
            value_size,
        }
    }

    /// `valueSize()`.
    pub fn value_size(&self) -> i32 {
        self.value_size
    }

    /// `intArray()`.
    pub fn int_array(&self) -> &dyn IntArray {
        self.int_array.as_ref()
    }

    /// `valueSize(int[] ia)` — min integer greater than all (non-negative) elements.
    pub fn value_size_of_slice(ia: &[i32]) -> i32 {
        let mut max = -1;
        for &i in ia {
            assert!(i >= 0, "{}", i);
            if i > max {
                max = i;
            }
        }
        max + 1
    }

    /// `valueSize(IntArray ia)`.
    pub fn value_size_of_arr(ia: &dyn IntArray) -> i32 {
        let mut max = -1;
        for j in 0..ia.size() {
            let i = ia.get(j);
            assert!(i >= 0, "{}", i);
            if i > max {
                max = i;
            }
        }
        max + 1
    }
}

impl IntArray for IndexArray {
    fn size(&self) -> i32 {
        self.int_array.size()
    }

    fn get(&self, index: i32) -> i32 {
        self.int_array.get(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::{create, packed_create};

    fn roundtrips(arr: &dyn IntArray, expect: &[i32]) {
        assert_eq!(arr.size(), expect.len() as i32);
        for (j, &v) in expect.iter().enumerate() {
            assert_eq!(arr.get(j as i32), v, "at index {j}");
        }
    }

    #[test]
    fn packed_backing_is_byte_exact() {
        // valueSize=4 -> 2 bits/value: [0,1,2,3] -> 0|1<<2|2<<4|3<<6 = 228
        let p = PackedIntArray::from_slice(&[0, 1, 2, 3], 4);
        assert_eq!(p.pack_index, 1);
        assert_eq!(p.ia, vec![228]);
        roundtrips(&p, &[0, 1, 2, 3]);

        // valueSize=2 -> 1 bit/value: eight 1s -> 0xFF
        let p = PackedIntArray::from_slice(&[1; 8], 2);
        assert_eq!(p.pack_index, 0);
        assert_eq!(p.ia, vec![255]);

        // valueSize=16 -> 4 bits/value
        let p = PackedIntArray::from_slice(&[15, 0, 0, 0], 16);
        assert_eq!(p.pack_index, 2);
        assert_eq!(p.ia, vec![15]);
    }

    #[test]
    fn packed_roundtrip_all_small_value_sizes() {
        for value_size in 1..=16 {
            let values: Vec<i32> = (0..40).map(|j| j % value_size).collect();
            let p = PackedIntArray::from_slice(&values, value_size);
            roundtrips(&p, &values);
            // IntList path must agree with the slice path.
            let mut il = IntList::new();
            values.iter().for_each(|&v| il.add(v));
            roundtrips(&PackedIntArray::from_list(&il, value_size), &values);
        }
    }

    #[test]
    fn packed_get_at_size_does_not_panic_ints1() {
        let p = PackedIntArray::from_slice(&[0, 1, 2, 3], 4);
        let _ = p.get(p.size()); // ints-1: `index > size`, so get(size) is allowed
    }

    #[test]
    #[should_panic]
    fn packed_get_past_size_panics() {
        let p = PackedIntArray::from_slice(&[0, 1, 2, 3], 4);
        let _ = p.get(p.size() + 1);
    }

    #[test]
    #[should_panic]
    fn packed_rejects_out_of_range_value() {
        let _ = PackedIntArray::from_slice(&[4], 4);
    }

    #[test]
    fn packed_byte_array_mask_quirk() {
        // signed path masks with 0xff; "unsigned" path masks with 0x7f (preserved).
        assert_eq!(
            PackedIntArray::from_signed_byte_array(&[0x83], 200).get(0),
            131
        );
        assert_eq!(
            PackedIntArray::from_unsigned_byte_array(&[0x83], 200).get(0),
            3
        );
    }

    #[test]
    fn packed_two_byte_array_big_endian() {
        let p = PackedIntArray::from_unsigned_two_byte_array(&[0x12, 0x34, 0x00, 0x05], 0x2000);
        roundtrips(&p, &[0x1234, 0x0005]);
    }

    #[test]
    fn unsigned_byte_array_basics() {
        roundtrips(
            &UnsignedByteArray::from_ints(&[0, 127, 255]),
            &[0, 127, 255],
        );
        roundtrips(&UnsignedByteArray::from_bytes(&[0, 1, 255]), &[0, 1, 255]);
        roundtrips(
            &UnsignedByteArray::from_ints_value_size(&[0, 5, 9], 10),
            &[0, 5, 9],
        );
    }

    #[test]
    #[should_panic]
    fn unsigned_byte_array_rejects_over_255() {
        let _ = UnsignedByteArray::from_ints(&[256]);
    }

    #[test]
    fn char_array_basics() {
        roundtrips(&CharArray::from_ints(&[0, 4660, 65535]), &[0, 4660, 65535]);
        // big-endian byte pairs: 0x1234 = 4660, 0x0005 = 5
        roundtrips(
            &CharArray::from_bytes(&[0x12, 0x34, 0x00, 0x05]),
            &[4660, 5],
        );
    }

    #[test]
    #[should_panic]
    fn char_array_odd_length_panics_ints2() {
        // ints-2: the `% 1` length check is dead, so odd input reads out of bounds.
        let _ = CharArray::from_bytes(&[0x12, 0x34, 0x00]);
    }

    #[test]
    fn wrapped_int_array_and_binary_search() {
        let w = WrappedIntArray::from_slice(&[1, 3, 5, 7, 9]);
        roundtrips(&w, &[1, 3, 5, 7, 9]);
        assert_eq!(w.binary_search(5), 2);
        assert_eq!(w.binary_search(4), -3); // insertion point 2 -> -(2)-1
        assert_eq!(w.binary_search(10), -6); // insertion point 5 -> -(5)-1
        roundtrips(&WrappedIntArray::from_int_stream(0..4), &[0, 1, 2, 3]);
    }

    #[test]
    fn index_array() {
        let ia = IndexArray::from_ints(&[0, 1, 2, 3], 4);
        roundtrips(&ia, &[0, 1, 2, 3]);
        assert_eq!(ia.value_size(), 4);
        assert_eq!(IndexArray::value_size_of_slice(&[3, 1, 2]), 4);
        assert_eq!(IndexArray::value_size_of_slice(&[]), 0);
    }

    #[test]
    #[should_panic]
    fn index_array_value_size_rejects_negative() {
        let _ = IndexArray::value_size_of_slice(&[-1]);
    }

    #[test]
    fn factories_route_and_roundtrip() {
        // packed_create: <=16 packed, <=256 byte, <=65536 char, else wrapped
        roundtrips(packed_create(&[0, 1, 2], 4).as_ref(), &[0, 1, 2]);
        roundtrips(packed_create(&[0, 200], 256).as_ref(), &[0, 200]);
        roundtrips(packed_create(&[0, 4660], 65536).as_ref(), &[0, 4660]);
        roundtrips(packed_create(&[0, 70000], 70001).as_ref(), &[0, 70000]);
        // create: 1/2/4 bytes only (no bit packing)
        roundtrips(create(&[0, 1], 2).as_ref(), &[0, 1]);
    }
}
