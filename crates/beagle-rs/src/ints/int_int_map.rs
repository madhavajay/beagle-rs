//! Port of `ints/IntIntMap.java` — a chained hash map with `int` keys and values.
//! Not thread-safe. Buckets hold a sorted (by key) singly-linked list of entries;
//! entries live in parallel `keys`/`values` arrays addressed via a `data` indirection.

use super::int_slice_to_string;

const NIL: i32 = -1;
const LOAD_FACTOR: f32 = 0.75;

/// `(int) Math.ceil(n / LOAD_FACTOR)`, reproducing Java's `float` division widened
/// to `double` for `Math.ceil`.
fn ceil_div_load(n: i32) -> i32 {
    ((n as f32 / LOAD_FACTOR) as f64).ceil() as i32
}

/// Port of `ints/IntIntMap.java`.
pub struct IntIntMap {
    size: i32,
    n_buckets: i32,
    next: Vec<i32>,
    data: Vec<i32>, // stores list index of keys and values
    keys: Vec<i32>,
    values: Vec<i32>,
    first_free_index: i32,
}

impl IntIntMap {
    /// `new IntIntMap(int capacity)`.
    pub fn new(capacity: i32) -> Self {
        assert!(capacity >= 1 && capacity <= (1 << 30), "{}", capacity);
        let num_buckets = ceil_div_load(capacity) + 1;
        let mut m = IntIntMap {
            size: 0,
            n_buckets: 0,
            next: Vec::new(),
            data: Vec::new(),
            keys: Vec::new(),
            values: Vec::new(),
            first_free_index: 0,
        };
        m.allocate_arrays(capacity, num_buckets);
        m.initialize_fields(num_buckets);
        m
    }

    fn allocate_arrays(&mut self, capacity: i32, num_buckets: i32) {
        let nb = num_buckets as usize;
        let cap = capacity as usize;
        self.next = vec![0; nb + cap];
        self.data = vec![0; nb + cap];
        self.keys = vec![0; cap];
        self.values = vec![0; cap];
    }

    fn initialize_fields(&mut self, num_buckets: i32) {
        self.size = 0;
        self.n_buckets = num_buckets;
        self.first_free_index = num_buckets;
        let nb = num_buckets as usize;
        for j in 0..nb {
            self.next[j] = NIL;
        }
        for j in nb..self.next.len() {
            self.next[j] = j as i32 + 1;
        }
    }

    fn rehash(&mut self, new_capacity: i32) {
        if new_capacity > self.size {
            let old_size = self.size;
            let old_keys = self.keys.clone();
            let old_values = self.values.clone();
            let new_num_buckets = ceil_div_load(new_capacity);
            self.allocate_arrays(new_capacity, new_num_buckets);
            self.initialize_fields(new_num_buckets);
            for j in 0..old_size as usize {
                self.put(old_keys[j], old_values[j]);
            }
        }
    }

    /// `clear()`.
    pub fn clear(&mut self) {
        let nb = self.n_buckets;
        self.initialize_fields(nb);
    }

    /// `contains(int key)`.
    pub fn contains(&self, key: i32) -> bool {
        self.index_of(key) >= 0
    }

    fn index_of(&self, key: i32) -> i32 {
        let mut index = self.next[self.bucket(key) as usize];
        while index != NIL && self.keys[self.data[index as usize] as usize] < key {
            index = self.next[index as usize];
        }
        if index != NIL && self.keys[self.data[index as usize] as usize] == key {
            index
        } else {
            -1
        }
    }

    /// `put(int key, int value)` — returns `true` if the map changed.
    pub fn put(&mut self, key: i32, value: i32) -> bool {
        let prev_index = self.prev_index(key);
        let next_index = self.next[prev_index as usize];
        if next_index == NIL || self.keys[self.data[next_index as usize] as usize] != key {
            let index = self.first_free_index;
            self.first_free_index = self.next[self.first_free_index as usize];
            self.next[prev_index as usize] = index;
            self.data[index as usize] = self.size;
            self.next[index as usize] = next_index;
            self.keys[self.size as usize] = key;
            self.values[self.size as usize] = value;
            self.size += 1;
            if self.size as usize == self.keys.len() {
                let new_capacity = 3 * self.keys.len() as i32 / 2 + 1;
                self.rehash(new_capacity);
            }
            true
        } else if self.values[self.data[next_index as usize] as usize] != value {
            let di = self.data[next_index as usize] as usize;
            self.values[di] = value;
            true
        } else {
            false
        }
    }

    /// `remove(int key)` — returns `true` if the map changed.
    pub fn remove(&mut self, key: i32) -> bool {
        let prev_index = self.prev_index(key);
        let index = self.next[prev_index as usize];
        if index == NIL || self.keys[self.data[index as usize] as usize] != key {
            false
        } else {
            let old_list_index = self.data[index as usize];
            self.next[prev_index as usize] = self.next[index as usize];
            self.next[index as usize] = self.first_free_index;
            self.first_free_index = index;

            self.size -= 1;
            if old_list_index != self.size {
                // overwrite removed key with the last list entry
                let moved_key = self.keys[self.size as usize];
                let idx = self.index_of(moved_key);
                self.data[idx as usize] = old_list_index;
                self.keys[old_list_index as usize] = self.keys[self.size as usize];
                self.values[old_list_index as usize] = self.values[self.size as usize];
            }
            true
        }
    }

    fn bucket(&self, key: i32) -> i32 {
        (71i32.wrapping_mul(key) % self.n_buckets).wrapping_abs()
    }

    fn prev_index(&self, key: i32) -> i32 {
        let mut prev_index = self.bucket(key);
        let mut index = self.next[prev_index as usize];
        while index != NIL && self.keys[self.data[index as usize] as usize] < key {
            prev_index = index;
            index = self.next[index as usize];
        }
        prev_index
    }

    /// `key(int index)`.
    pub fn key(&self, index: i32) -> i32 {
        assert!(index < self.size, "{}", index);
        self.keys[index as usize]
    }

    /// `get(int key, int sentinel)`.
    pub fn get(&self, key: i32, sentinel: i32) -> i32 {
        let index = self.index_of(key);
        if index == -1 {
            sentinel
        } else {
            self.values[self.data[index as usize] as usize]
        }
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.size
    }

    /// `keys()` — the keys in current index order.
    pub fn keys(&self) -> Vec<i32> {
        self.keys[0..self.size as usize].to_vec()
    }

    /// `values()` — the values in current key-index order.
    pub fn values(&self) -> Vec<i32> {
        self.values[0..self.size as usize].to_vec()
    }
}

impl std::fmt::Display for IntIntMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&int_slice_to_string(&self.keys()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn put_get_contains_remove() {
        let mut m = IntIntMap::new(4);
        assert!(m.put(10, 100));
        assert!(m.put(20, 200));
        assert!(!m.put(10, 100)); // unchanged
        assert!(m.put(10, 101)); // value changed
        assert_eq!(m.get(10, -1), 101);
        assert_eq!(m.get(20, -1), 200);
        assert_eq!(m.get(30, -1), -1); // sentinel
        assert!(m.contains(10));
        assert!(!m.contains(30));
        assert_eq!(m.size(), 2);
        assert!(m.remove(10));
        assert!(!m.remove(10));
        assert_eq!(m.get(10, -1), -1);
        assert_eq!(m.size(), 1);
        m.clear();
        assert_eq!(m.size(), 0);
    }

    #[test]
    fn rehash_preserves_entries() {
        let mut m = IntIntMap::new(4);
        for k in 0..100 {
            m.put(k, k * 3);
        }
        assert_eq!(m.size(), 100);
        for k in 0..100 {
            assert_eq!(m.get(k, -1), k * 3);
        }
    }

    #[test]
    #[should_panic]
    fn rejects_zero_capacity_ints7() {
        let _ = IntIntMap::new(0);
    }

    #[test]
    fn matches_hashmap_over_deterministic_ops() {
        let mut m = IntIntMap::new(4);
        let mut h: HashMap<i32, i32> = HashMap::new();
        for i in 0..5000i32 {
            let key = i % 100;
            let value = i.wrapping_mul(2_654_435_761u32 as i32);
            match i % 200 {
                0 => {
                    m.clear();
                    h.clear();
                }
                r if r % 3 == 1 => {
                    m.put(key, value);
                    h.insert(key, value);
                }
                r if r % 3 == 2 => {
                    m.remove(key);
                    h.remove(&key);
                }
                _ => {
                    let got = m.get(key, -1);
                    match h.get(&key) {
                        Some(&v) => assert_eq!(got, v),
                        None => assert_eq!(got, -1),
                    }
                }
            }
            assert_eq!(m.size() as usize, h.len(), "size mismatch at i={i}");
        }
    }
}
