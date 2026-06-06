//! Port of `ints/IntList.java` and `ints/SynchedIntList.java`.

use super::{int_slice_to_string, java_binary_search};

/// Port of `ints/IntList.java` — a growable list of `int`, with `clear()` but no
/// `remove()`. Not thread-safe. The backing store grows like the Java version and
/// is never cleared on `pop`/`truncate`/`clear`, so `copy_of`/`copy_of_range` see
/// the same (occasionally stale) backing slots the Java code would.
#[derive(Clone)]
pub struct IntList {
    size: i32,
    values: Vec<i32>,
}

/// `IntList.DEFAULT_INIT_CAPACITY`.
pub const DEFAULT_INIT_CAPACITY: i32 = 16;

impl IntList {
    /// `new IntList()`.
    pub fn new() -> Self {
        Self::with_capacity(DEFAULT_INIT_CAPACITY)
    }

    /// `new IntList(int initCapacity)`.
    pub fn with_capacity(init_capacity: i32) -> Self {
        assert!(init_capacity >= 0, "{}", init_capacity);
        IntList {
            size: 0,
            values: vec![0; init_capacity as usize],
        }
    }

    /// `new IntList(int[] ia)` — clones the array; capacity equals length.
    pub fn from_slice(ia: &[i32]) -> Self {
        IntList {
            size: ia.len() as i32,
            values: ia.to_vec(),
        }
    }

    /// `add(int value)`.
    pub fn add(&mut self, value: i32) {
        if self.size as usize == self.values.len() {
            let new_capacity = (self.values.len() * 3) / 2 + 1;
            self.values.resize(new_capacity, 0);
        }
        self.values[self.size as usize] = value;
        self.size += 1;
    }

    /// `pop()` — removes and returns the last entry (does not clear the slot).
    pub fn pop(&mut self) -> i32 {
        self.size -= 1;
        self.values[self.size as usize]
    }

    /// `get(int index)` — throws (panics) only when `index >= size` (ints-5).
    pub fn get(&self, index: i32) -> i32 {
        assert!(index < self.size, "{}", index);
        self.values[index as usize]
    }

    /// `set(int index, int value)` — returns the previous value.
    pub fn set(&mut self, index: i32, value: i32) -> i32 {
        assert!(index < self.size, "{}", index);
        let old = self.values[index as usize];
        self.values[index as usize] = value;
        old
    }

    /// `getAndIncrement(int index)`.
    pub fn get_and_increment(&mut self, index: i32) -> i32 {
        assert!(index >= 0 && index < self.size, "{}", index);
        let old = self.values[index as usize];
        self.values[index as usize] = old.wrapping_add(1);
        old
    }

    /// `getAndDecrement(int index)`.
    pub fn get_and_decrement(&mut self, index: i32) -> i32 {
        assert!(index >= 0 && index < self.size, "{}", index);
        let old = self.values[index as usize];
        self.values[index as usize] = old.wrapping_sub(1);
        old
    }

    /// `incrementAndGet(int index)`.
    pub fn increment_and_get(&mut self, index: i32) -> i32 {
        assert!(index >= 0 && index < self.size, "{}", index);
        let v = self.values[index as usize].wrapping_add(1);
        self.values[index as usize] = v;
        v
    }

    /// `decrementAndGet(int index)`.
    pub fn decrement_and_get(&mut self, index: i32) -> i32 {
        assert!(index >= 0 && index < self.size, "{}", index);
        let v = self.values[index as usize].wrapping_sub(1);
        self.values[index as usize] = v;
        v
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.size
    }

    /// `isEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// `sort()` — ascending.
    pub fn sort(&mut self) {
        let n = self.size as usize;
        self.values[0..n].sort_unstable();
    }

    /// `binarySearch(int value)` — Java semantics over `[0, size)`.
    pub fn binary_search(&self, value: i32) -> i32 {
        java_binary_search(&self.values, 0, self.size, value)
    }

    /// `copyOf(int newLength)` — truncate/pad-with-0 copy of the logical contents.
    pub fn copy_of(&self, new_length: i32) -> Vec<i32> {
        let mut out = vec![0i32; new_length as usize];
        let copy_n = new_length.min(self.size).max(0) as usize;
        out[0..copy_n].copy_from_slice(&self.values[0..copy_n]);
        out
    }

    /// `copyOfRange(int start, int end)` — see ints-4 for the negative-length quirk.
    pub fn copy_of_range(&self, start: i32, end: i32) -> Vec<i32> {
        assert!(start >= 0 && start <= self.size, "{}", start);
        let len = end - start; // NegativeArraySizeException analog if end < start
        let mut out = vec![0i32; len as usize];
        let n = if end < self.size {
            (end - start) as usize
        } else {
            (self.size - start) as usize
        };
        out[0..n].copy_from_slice(&self.values[start as usize..start as usize + n]);
        out
    }

    /// `truncate(int newSize)`.
    pub fn truncate(&mut self, new_size: i32) {
        assert!(new_size >= 0, "{}", new_size);
        if new_size < self.size {
            self.size = new_size;
        }
    }

    /// `toArray()` — the logical contents.
    pub fn to_array(&self) -> Vec<i32> {
        self.values[0..self.size as usize].to_vec()
    }

    /// `clear()`.
    pub fn clear(&mut self) {
        self.size = 0;
    }
}

impl Default for IntList {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Display for IntList {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&int_slice_to_string(&self.to_array()))
    }
}

/// Port of `ints/SynchedIntList.java` — a thread-safe `IntList` subset. Implemented
/// by guarding an [`IntList`] with a mutex so all methods take `&self`, mirroring
/// Java's `synchronized` instance methods.
pub struct SynchedIntList {
    inner: std::sync::Mutex<IntList>,
}

impl SynchedIntList {
    /// `new SynchedIntList()`.
    pub fn new() -> Self {
        Self::with_capacity(DEFAULT_INIT_CAPACITY)
    }

    /// `new SynchedIntList(int initCapacity)`.
    pub fn with_capacity(init_capacity: i32) -> Self {
        SynchedIntList {
            inner: std::sync::Mutex::new(IntList::with_capacity(init_capacity)),
        }
    }

    /// `new SynchedIntList(int[] ia)`.
    pub fn from_slice(ia: &[i32]) -> Self {
        SynchedIntList {
            inner: std::sync::Mutex::new(IntList::from_slice(ia)),
        }
    }

    /// `add(int value)`.
    pub fn add(&self, value: i32) {
        self.inner.lock().unwrap().add(value);
    }

    /// `get(int index)`.
    pub fn get(&self, index: i32) -> i32 {
        self.inner.lock().unwrap().get(index)
    }

    /// `set(int index, int value)`.
    pub fn set(&self, index: i32, value: i32) -> i32 {
        self.inner.lock().unwrap().set(index, value)
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.inner.lock().unwrap().size()
    }

    /// `toArray()`.
    pub fn to_array(&self) -> Vec<i32> {
        self.inner.lock().unwrap().to_array()
    }

    /// `clear()`.
    pub fn clear(&self) {
        self.inner.lock().unwrap().clear();
    }
}

impl Default for SynchedIntList {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Display for SynchedIntList {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&int_slice_to_string(&self.to_array()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_get_set_grow() {
        let mut l = IntList::new();
        for j in 0..20 {
            l.add(j * 10);
        }
        assert_eq!(l.size(), 20);
        assert_eq!(l.get(0), 0);
        assert_eq!(l.get(19), 190);
        assert_eq!(l.set(5, 999), 50);
        assert_eq!(l.get(5), 999);
    }

    #[test]
    fn pop_does_not_clear_backing() {
        let mut l = IntList::from_slice(&[1, 2, 3]);
        assert_eq!(l.pop(), 3);
        assert_eq!(l.size(), 2);
        l.add(7);
        assert_eq!(l.to_array(), vec![1, 2, 7]);
    }

    #[test]
    fn copy_of_truncates_and_pads() {
        let l = IntList::from_slice(&[1, 2, 3]);
        assert_eq!(l.copy_of(2), vec![1, 2]);
        assert_eq!(l.copy_of(5), vec![1, 2, 3, 0, 0]);
    }

    #[test]
    fn copy_of_range_within_and_past_size() {
        let l = IntList::from_slice(&[10, 11, 12, 13, 14]);
        assert_eq!(l.copy_of_range(1, 3), vec![11, 12]);
        assert_eq!(l.copy_of_range(3, 7), vec![13, 14, 0, 0]); // end > size pads
    }

    #[test]
    #[should_panic]
    fn copy_of_range_negative_length_panics_ints4() {
        let l = IntList::from_slice(&[10, 11, 12]);
        let _ = l.copy_of_range(2, 1);
    }

    #[test]
    fn truncate_sort_binary_search() {
        let mut l = IntList::from_slice(&[5, 1, 4, 2, 3]);
        l.sort();
        assert_eq!(l.to_array(), vec![1, 2, 3, 4, 5]);
        assert_eq!(l.binary_search(4), 3);
        assert_eq!(l.binary_search(0), -1); // insertion point 0 -> -(0)-1
        l.truncate(2);
        assert_eq!(l.to_array(), vec![1, 2]);
        l.truncate(10); // no-op when newSize > size
        assert_eq!(l.size(), 2);
    }

    #[test]
    fn increment_decrement_family() {
        let mut l = IntList::from_slice(&[10]);
        assert_eq!(l.get_and_increment(0), 10);
        assert_eq!(l.get(0), 11);
        assert_eq!(l.increment_and_get(0), 12);
        assert_eq!(l.get_and_decrement(0), 12);
        assert_eq!(l.decrement_and_get(0), 10);
    }

    #[test]
    fn display_and_clear() {
        let mut l = IntList::from_slice(&[1, 2, 3]);
        assert_eq!(l.to_string(), "[1, 2, 3]");
        l.clear();
        assert_eq!(l.size(), 0);
        assert_eq!(IntList::new().to_string(), "[]");
    }

    #[test]
    #[should_panic]
    fn get_out_of_bounds_panics() {
        let l = IntList::from_slice(&[1, 2, 3]);
        let _ = l.get(3);
    }

    #[test]
    fn synched_int_list_basics() {
        let l = SynchedIntList::new();
        l.add(1);
        l.add(2);
        assert_eq!(l.size(), 2);
        assert_eq!(l.set(0, 9), 1);
        assert_eq!(l.get(0), 9);
        assert_eq!(l.to_array(), vec![9, 2]);
        l.clear();
        assert_eq!(l.size(), 0);
    }
}
