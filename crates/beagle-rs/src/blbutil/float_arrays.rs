//! Port of `blbutil/FloatArray.java`, `blbutil/DoubleArray.java`, and
//! `blbutil/FloatList.java`.
//!
//! Note: the Java `toString()` of these types uses `Float`/`Double.toString`, whose
//! exact formatting differs from Rust's and is documented as "unspecified and subject
//! to change". The `Display` impls here are best-effort and are **not** on the
//! byte-for-byte output path (VCF float formatting lives in `vcf::VcfWriter`).

use super::{binary_search_f32, binary_search_f64};

/// Port of `blbutil/FloatArray.java` — an immutable list of `float`.
#[derive(Clone)]
pub struct FloatArray {
    values: Vec<f32>,
}

impl FloatArray {
    /// `new FloatArray(float[] values)`.
    pub fn from_floats(values: &[f32]) -> Self {
        FloatArray {
            values: values.to_vec(),
        }
    }

    /// `new FloatArray(double[] values)` — narrows each element to `float`.
    pub fn from_doubles(values: &[f64]) -> Self {
        FloatArray {
            values: values.iter().map(|&d| d as f32).collect(),
        }
    }

    /// `FloatArray.fromIntBits(int[] bits)` — `Float.intBitsToFloat` of each element.
    pub fn from_int_bits(bits: &[i32]) -> Self {
        FloatArray {
            values: bits.iter().map(|&b| f32::from_bits(b as u32)).collect(),
        }
    }

    /// `get(int index)`.
    pub fn get(&self, index: i32) -> f32 {
        self.values[index as usize]
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.values.len() as i32
    }

    /// `isEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// `binarySearch(float key)`.
    pub fn binary_search(&self, key: f32) -> i32 {
        binary_search_f32(&self.values, 0, self.values.len() as i32, key)
    }

    /// `binarySearch(int fromIndex, int toIndex, float key)`.
    pub fn binary_search_range(&self, from_index: i32, to_index: i32, key: f32) -> i32 {
        assert!(
            from_index <= to_index,
            "fromIndex({from_index}) > toIndex({to_index})"
        );
        binary_search_f32(&self.values, from_index, to_index, key)
    }

    /// `toArray()`.
    pub fn to_array(&self) -> Vec<f32> {
        self.values.clone()
    }
}

/// Port of `blbutil/DoubleArray.java` — an immutable list of `double`.
#[derive(Clone)]
pub struct DoubleArray {
    values: Vec<f64>,
}

impl DoubleArray {
    /// `new DoubleArray(double[] values)`.
    pub fn from_doubles(values: &[f64]) -> Self {
        DoubleArray {
            values: values.to_vec(),
        }
    }

    /// `new DoubleArray(DoubleStream values)`.
    pub fn from_double_stream<I: IntoIterator<Item = f64>>(values: I) -> Self {
        DoubleArray {
            values: values.into_iter().collect(),
        }
    }

    /// `get(int index)`.
    pub fn get(&self, index: i32) -> f64 {
        self.values[index as usize]
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.values.len() as i32
    }

    /// `isEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// `binarySearch(double key)`.
    pub fn binary_search(&self, key: f64) -> i32 {
        binary_search_f64(&self.values, 0, self.values.len() as i32, key)
    }

    /// `binarySearch(int fromIndex, int toIndex, double key)`.
    pub fn binary_search_range(&self, from_index: i32, to_index: i32, key: f64) -> i32 {
        assert!(
            from_index <= to_index,
            "fromIndex({from_index}) > toIndex({to_index})"
        );
        binary_search_f64(&self.values, from_index, to_index, key)
    }

    /// `toArray()`.
    pub fn to_array(&self) -> Vec<f64> {
        self.values.clone()
    }
}

/// Port of `blbutil/FloatList.java` — a growable list of `float` (`clear` but no
/// `remove`). Default initial capacity 10.
pub struct FloatList {
    size: i32,
    values: Vec<f32>,
}

/// `FloatList.DEFAULT_INIT_CAPACITY`.
pub const FLOAT_LIST_DEFAULT_INIT_CAPACITY: i32 = 10;

impl FloatList {
    /// `new FloatList()`.
    pub fn new() -> Self {
        Self::with_capacity(FLOAT_LIST_DEFAULT_INIT_CAPACITY)
    }

    /// `new FloatList(int initCapacity)`.
    pub fn with_capacity(init_capacity: i32) -> Self {
        assert!(init_capacity >= 0, "{}", init_capacity);
        FloatList {
            size: 0,
            values: vec![0.0; init_capacity as usize],
        }
    }

    /// `add(float element)`.
    pub fn add(&mut self, element: f32) {
        if self.size as usize == self.values.len() {
            let new_capacity = (self.values.len() * 3) / 2 + 1;
            self.values.resize(new_capacity, 0.0);
        }
        self.values[self.size as usize] = element;
        self.size += 1;
    }

    /// `addToElement(int index, float value)`.
    pub fn add_to_element(&mut self, index: i32, value: f32) {
        assert!(index < self.size, "{}", index);
        self.values[index as usize] += value;
    }

    /// `get(int index)`.
    pub fn get(&self, index: i32) -> f32 {
        assert!(index < self.size, "{}", index);
        self.values[index as usize]
    }

    /// `set(int index, float value)` — returns the previous value.
    pub fn set(&mut self, index: i32, value: f32) -> f32 {
        assert!(index < self.size, "{}", index);
        let old = self.values[index as usize];
        self.values[index as usize] = value;
        old
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.size
    }

    /// `isEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// `toArray()`.
    pub fn to_array(&self) -> Vec<f32> {
        self.values[0..self.size as usize].to_vec()
    }

    /// `clear()`.
    pub fn clear(&mut self) {
        self.size = 0;
    }
}

impl Default for FloatList {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_array_roundtrip_and_search() {
        let a = FloatArray::from_floats(&[1.5, 2.5, 3.5]);
        assert_eq!(a.size(), 3);
        assert_eq!(a.get(1), 2.5);
        assert_eq!(a.binary_search(3.5), 2);
        assert_eq!(a.to_array(), vec![1.5, 2.5, 3.5]);
        assert!(!a.is_empty());
    }

    #[test]
    fn float_array_from_int_bits() {
        let bits = [1.0f32.to_bits() as i32, 2.0f32.to_bits() as i32];
        let a = FloatArray::from_int_bits(&bits);
        assert_eq!(a.get(0), 1.0);
        assert_eq!(a.get(1), 2.0);
    }

    #[test]
    fn double_array_basics() {
        let a = DoubleArray::from_doubles(&[1.0, 2.0, 4.0]);
        assert_eq!(a.size(), 3);
        assert_eq!(a.get(2), 4.0);
        assert_eq!(a.binary_search(2.0), 1);
        let b = DoubleArray::from_double_stream([5.0, 6.0]);
        assert_eq!(b.to_array(), vec![5.0, 6.0]);
    }

    #[test]
    fn float_list_grow_and_ops() {
        let mut l = FloatList::new();
        for j in 0..15 {
            l.add(j as f32);
        }
        assert_eq!(l.size(), 15);
        assert_eq!(l.get(14), 14.0);
        assert_eq!(l.set(0, 99.0), 0.0);
        l.add_to_element(0, 1.0);
        assert_eq!(l.get(0), 100.0);
        l.clear();
        assert!(l.is_empty());
    }

    #[test]
    #[should_panic]
    fn float_list_get_oob_panics() {
        let l = FloatList::new();
        let _ = l.get(0);
    }
}
