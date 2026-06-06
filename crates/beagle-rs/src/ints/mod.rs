//! Port of the Java `ints` package: immutable integer-packed arrays and an
//! int→int map. Mirrors `ints/IntArray.java` (this module) plus one Rust module
//! per implementing class.
//!
//! Faithfulness notes:
//! - Indices and sizes are `i32` to match Java `int` semantics exactly (including
//!   the documented out-of-bounds quirks; see `docs/known-quirks.md`).
//! - Validation failures that throw `IllegalArgumentException` / bounds exceptions
//!   in Java are reproduced as `panic!` here (these guard internal invariants).

mod arrays;
mod int_int_map;
mod lists;

pub use arrays::{CharArray, IndexArray, PackedIntArray, UnsignedByteArray, WrappedIntArray};
pub use int_int_map::IntIntMap;
pub use lists::{IntList, SynchedIntList};

/// Port of the Java `IntArray` interface — an immutable `int[]`.
pub trait IntArray {
    /// Returns the number of elements.
    fn size(&self) -> i32;

    /// Returns the element at `index`. Out-of-bounds behavior matches the Java
    /// implementor exactly (some implementors have off-by-one quirks).
    fn get(&self, index: i32) -> i32;
}

// ---- Java `int` shift semantics (the JVM masks the shift count to 5 bits) ----

/// Java `x << n` on `int`.
#[inline]
pub(crate) fn jshl(x: i32, n: i32) -> i32 {
    x.wrapping_shl(n as u32)
}

/// Java `x >> n` on `int` (arithmetic / sign-propagating).
#[inline]
pub(crate) fn jshr(x: i32, n: i32) -> i32 {
    x.wrapping_shr(n as u32)
}

/// Java `x >>> n` on `int` (logical / zero-fill).
#[inline]
pub(crate) fn jushr(x: i32, n: i32) -> i32 {
    (x as u32).wrapping_shr(n as u32) as i32
}

/// Java `java.util.Arrays.binarySearch(a, fromIndex, toIndex, key)` on a sorted
/// ascending range. Returns the index if found, else `-(insertionPoint) - 1`.
pub(crate) fn java_binary_search(a: &[i32], from: i32, to: i32, key: i32) -> i32 {
    let mut low = from;
    let mut high = to - 1;
    while low <= high {
        let mid = ((low as u32).wrapping_add(high as u32) >> 1) as i32; // (low+high) >>> 1
        let midval = a[mid as usize];
        if midval < key {
            low = mid + 1;
        } else if midval > key {
            high = mid - 1;
        } else {
            return mid;
        }
    }
    -(low + 1)
}

// ---- static helpers from IntArray.java ----

/// `IntArray.toArray` — a copy of `ia` as a `Vec<i32>`.
pub fn to_array(ia: &dyn IntArray) -> Vec<i32> {
    let n = ia.size();
    let mut copy = Vec::with_capacity(n.max(0) as usize);
    for j in 0..n {
        copy.push(ia.get(j));
    }
    copy
}

/// `IntArray.asString` — `java.util.Arrays.toString` of the equivalent `int[]`.
pub fn as_string(ia: &dyn IntArray) -> String {
    int_slice_to_string(&to_array(ia))
}

/// `java.util.Arrays.toString(int[])`: `[]`, `[1]`, `[1, 2, 3]`.
pub(crate) fn int_slice_to_string(a: &[i32]) -> String {
    let mut s = String::from("[");
    for (i, v) in a.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        s.push_str(&v.to_string());
    }
    s.push(']');
    s
}

/// `IntArray.equals` — same length and element sequence.
pub fn equals(a: &dyn IntArray, b: &dyn IntArray) -> bool {
    if a.size() != b.size() {
        return false;
    }
    for j in 0..a.size() {
        if a.get(j) != b.get(j) {
            return false;
        }
    }
    true
}

/// `IntArray.max` — maximum element, or `i32::MIN` if empty.
pub fn max(ia: &dyn IntArray) -> i32 {
    let mut m = i32::MIN;
    for j in 0..ia.size() {
        let v = ia.get(j);
        if v > m {
            m = v;
        }
    }
    m
}

/// `IntArray.min` — minimum element, or `i32::MAX` if empty.
pub fn min(ia: &dyn IntArray) -> i32 {
    let mut m = i32::MAX;
    for j in 0..ia.size() {
        let v = ia.get(j);
        if v < m {
            m = v;
        }
    }
    m
}

/// `IntArray.packedCreate(int[], int)` — packs into 1/2/4 bits, then 1/2/4 bytes.
pub fn packed_create(ia: &[i32], value_size: i32) -> Box<dyn IntArray> {
    assert!(value_size >= 1, "{}", value_size);
    if value_size <= 16 {
        Box::new(PackedIntArray::from_slice(ia, value_size))
    } else if value_size <= 256 {
        Box::new(UnsignedByteArray::from_ints_value_size(ia, value_size))
    } else if value_size <= 65536 {
        Box::new(CharArray::from_ints_value_size(ia, value_size))
    } else {
        Box::new(WrappedIntArray::from_slice_value_size(ia, value_size))
    }
}

/// `IntArray.packedCreate(IntList, int)`.
pub fn packed_create_from_list(il: &IntList, value_size: i32) -> Box<dyn IntArray> {
    assert!(value_size >= 1, "{}", value_size);
    if value_size <= 16 {
        Box::new(PackedIntArray::from_list(il, value_size))
    } else if value_size <= 256 {
        Box::new(UnsignedByteArray::from_list_value_size(il, value_size))
    } else if value_size <= 65536 {
        Box::new(CharArray::from_list_value_size(il, value_size))
    } else {
        Box::new(WrappedIntArray::from_list_value_size(il, value_size))
    }
}

/// `IntArray.create(int[], int)` — stored in 1, 2, or 4 bytes.
pub fn create(ia: &[i32], value_size: i32) -> Box<dyn IntArray> {
    assert!(value_size >= 1, "{}", value_size);
    if value_size <= 256 {
        Box::new(UnsignedByteArray::from_ints_value_size(ia, value_size))
    } else if value_size <= 65536 {
        Box::new(CharArray::from_ints_value_size(ia, value_size))
    } else {
        Box::new(WrappedIntArray::from_slice_value_size(ia, value_size))
    }
}

/// `IntArray.create(IntList, int)`.
pub fn create_from_list(il: &IntList, value_size: i32) -> Box<dyn IntArray> {
    assert!(value_size >= 1, "{}", value_size);
    if value_size <= 256 {
        Box::new(UnsignedByteArray::from_list_value_size(il, value_size))
    } else if value_size <= 65536 {
        Box::new(CharArray::from_list_value_size(il, value_size))
    } else {
        Box::new(WrappedIntArray::from_list_value_size(il, value_size))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn statics_over_intarray() {
        let a = WrappedIntArray::from_slice(&[3, 1, 4, 1, 5]);
        assert_eq!(to_array(&a), vec![3, 1, 4, 1, 5]);
        assert_eq!(as_string(&a), "[3, 1, 4, 1, 5]");
        assert_eq!(max(&a), 5);
        assert_eq!(min(&a), 1);

        let b = WrappedIntArray::from_slice(&[3, 1, 4, 1, 5]);
        let c = WrappedIntArray::from_slice(&[3, 1, 4, 1, 6]);
        assert!(equals(&a, &b));
        assert!(!equals(&a, &c));
    }

    #[test]
    fn statics_on_empty() {
        let e = WrappedIntArray::from_slice(&[]);
        assert_eq!(as_string(&e), "[]");
        assert_eq!(max(&e), i32::MIN);
        assert_eq!(min(&e), i32::MAX);
    }

    #[test]
    fn binary_search_semantics() {
        let a = [1, 3, 5, 7, 9];
        assert_eq!(java_binary_search(&a, 0, 5, 7), 3);
        assert_eq!(java_binary_search(&a, 0, 5, 2), -2); // insertion point 1
        assert_eq!(java_binary_search(&a, 0, 5, 0), -1);
        assert_eq!(java_binary_search(&a, 0, 5, 10), -6);
    }
}
