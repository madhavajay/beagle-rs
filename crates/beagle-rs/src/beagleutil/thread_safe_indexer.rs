//! Port of `beagleutil/ThreadSafeIndexer.java` — assigns consecutive indices to
//! distinct objects. Java synchronizes every method; here a `Mutex` guards the state.

use std::collections::HashMap;
use std::hash::Hash;
use std::sync::Mutex;

/// `ThreadSafeIndexer.DEFAULT_INIT_CAPACITY`.
pub const DEFAULT_INIT_CAPACITY: i32 = 500;

struct Inner<T> {
    list: Vec<T>,
    map: HashMap<T, i32>,
}

/// Port of `beagleutil/ThreadSafeIndexer.java`.
pub struct ThreadSafeIndexer<T: Eq + Hash + Clone> {
    inner: Mutex<Inner<T>>,
}

impl<T: Eq + Hash + Clone> ThreadSafeIndexer<T> {
    /// `new ThreadSafeIndexer()`.
    pub fn new() -> Self {
        Self::with_capacity(DEFAULT_INIT_CAPACITY)
    }

    /// `new ThreadSafeIndexer(int initCapacity)`.
    pub fn with_capacity(init_capacity: i32) -> Self {
        assert!(init_capacity >= 1, "{}", init_capacity);
        ThreadSafeIndexer {
            inner: Mutex::new(Inner {
                list: Vec::with_capacity(init_capacity as usize),
                map: HashMap::with_capacity(init_capacity as usize),
            }),
        }
    }

    /// `getIndex(T object)` — index the object if needed, return its index.
    pub fn get_index(&self, object: T) -> i32 {
        let mut g = self.inner.lock().unwrap();
        if let Some(&i) = g.map.get(&object) {
            return i;
        }
        let idx = g.list.len() as i32;
        g.list.push(object.clone());
        g.map.insert(object, idx);
        idx
    }

    /// `getIndices(T[] objects)`.
    pub fn get_indices(&self, objects: &[T]) -> Vec<i32> {
        let mut g = self.inner.lock().unwrap();
        let mut indices = Vec::with_capacity(objects.len());
        for object in objects {
            if let Some(&i) = g.map.get(object) {
                indices.push(i);
            } else {
                let idx = g.list.len() as i32;
                g.list.push(object.clone());
                g.map.insert(object.clone(), idx);
                indices.push(idx);
            }
        }
        indices
    }

    /// `getIndexIfIndexed(T object)` — index, or `-1` if not present.
    pub fn get_index_if_indexed(&self, object: &T) -> i32 {
        let g = self.inner.lock().unwrap();
        g.map.get(object).copied().unwrap_or(-1)
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.inner.lock().unwrap().list.len() as i32
    }

    /// `item(int index)`.
    pub fn item(&self, index: i32) -> T {
        self.inner.lock().unwrap().list[index as usize].clone()
    }

    /// `items(int[] indices)`.
    pub fn items_at(&self, indices: &[i32]) -> Vec<T> {
        let g = self.inner.lock().unwrap();
        indices
            .iter()
            .map(|&i| g.list[i as usize].clone())
            .collect()
    }

    /// `items()`.
    pub fn items(&self) -> Vec<T> {
        self.inner.lock().unwrap().list.clone()
    }
}

impl<T: Eq + Hash + Clone> Default for ThreadSafeIndexer<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assigns_consecutive_indices() {
        let ix: ThreadSafeIndexer<String> = ThreadSafeIndexer::new();
        assert_eq!(ix.get_index("a".into()), 0);
        assert_eq!(ix.get_index("b".into()), 1);
        assert_eq!(ix.get_index("a".into()), 0); // existing -> same index
        assert_eq!(ix.get_index("c".into()), 2);
        assert_eq!(ix.size(), 3);
        assert_eq!(ix.item(1), "b");
        assert_eq!(ix.get_index_if_indexed(&"c".into()), 2);
        assert_eq!(ix.get_index_if_indexed(&"z".into()), -1);
        assert_eq!(
            ix.get_indices(&["b".into(), "d".into(), "a".into()]),
            vec![1, 3, 0]
        );
        assert_eq!(ix.items_at(&[3, 0]), vec!["d".to_string(), "a".to_string()]);
        assert_eq!(ix.items(), vec!["a", "b", "c", "d"]);
    }

    #[test]
    #[should_panic]
    fn rejects_capacity_below_one() {
        let _: ThreadSafeIndexer<String> = ThreadSafeIndexer::with_capacity(0);
    }
}
