//! Faithful reimplementations of JDK standard-library classes that Beagle uses
//! directly (not part of any Beagle package). Determinism of the port depends on
//! these matching the JDK bit-for-bit.

use std::cmp::Ordering;

/// Port of `java.util.Random` — the 48-bit linear congruential generator. Beagle
/// seeds it deterministically from `seed=`, so reproducing it exactly is required for
/// byte-for-byte phasing/imputation parity.
pub struct Random {
    seed: i64,
}

const MULTIPLIER: i64 = 0x5DEECE66D;
const ADDEND: i64 = 0xB;
const MASK: i64 = (1 << 48) - 1;
const DOUBLE_UNIT: f64 = 1.0 / ((1i64 << 53) as f64); // 0x1.0p-53
const FLOAT_UNIT: f32 = 1.0 / ((1i32 << 24) as f32); // 0x1.0p-24

impl Random {
    /// `new Random(long seed)`.
    pub fn new(seed: i64) -> Self {
        Random {
            seed: (seed ^ MULTIPLIER) & MASK,
        }
    }

    /// `setSeed(long seed)`.
    pub fn set_seed(&mut self, seed: i64) {
        self.seed = (seed ^ MULTIPLIER) & MASK;
    }

    /// `next(int bits)` — the protected core generator.
    fn next(&mut self, bits: i32) -> i32 {
        self.seed = self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND) & MASK;
        // (int)(seed >>> (48 - bits)); seed is non-negative (masked), so >> is logical.
        (self.seed >> (48 - bits)) as i32
    }

    /// `nextInt()`.
    pub fn next_int(&mut self) -> i32 {
        self.next(32)
    }

    /// `nextInt(int bound)`.
    pub fn next_int_bound(&mut self, bound: i32) -> i32 {
        assert!(bound > 0, "bound must be positive: {}", bound);
        if (bound & bound.wrapping_neg()) == bound {
            // bound is a power of two
            return ((bound as i64).wrapping_mul(self.next(31) as i64) >> 31) as i32;
        }
        loop {
            let bits = self.next(31);
            let val = bits % bound;
            if bits.wrapping_sub(val).wrapping_add(bound - 1) >= 0 {
                return val;
            }
        }
    }

    /// `nextLong()`.
    pub fn next_long(&mut self) -> i64 {
        let hi = (self.next(32) as i64) << 32;
        hi.wrapping_add(self.next(32) as i64)
    }

    /// `nextBoolean()`.
    pub fn next_boolean(&mut self) -> bool {
        self.next(1) != 0
    }

    /// `nextDouble()`.
    pub fn next_double(&mut self) -> f64 {
        let hi = (self.next(26) as i64) << 27;
        (hi.wrapping_add(self.next(27) as i64)) as f64 * DOUBLE_UNIT
    }

    /// `nextFloat()`.
    pub fn next_float(&mut self) -> f32 {
        self.next(24) as f32 * FLOAT_UNIT
    }
}

/// Port of `java.util.PriorityQueue` — an array-backed binary min-heap ordered by `Ord`
/// (mirroring elements' `compareTo`). The exact `siftUp`/`siftDown` element movement is
/// reproduced so that ties (elements that compare `Equal`) are broken identically to the JDK,
/// which Beagle's composite-haplotype construction depends on for byte-for-byte output.
///
/// Only the operations Beagle uses are provided: `offer`/`add`, `poll`, `peek`, `clear`,
/// `size`, `is_empty`. (`remove(Object)`/`contains` are unused and omitted.)
pub struct PriorityQueue<T: Ord> {
    queue: Vec<T>,
}

impl<T: Ord> PriorityQueue<T> {
    /// `new PriorityQueue(int initialCapacity)`.
    pub fn new(initial_capacity: usize) -> Self {
        PriorityQueue {
            queue: Vec::with_capacity(initial_capacity),
        }
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.queue.len() as i32
    }

    /// `isEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }

    /// `clear()`.
    pub fn clear(&mut self) {
        self.queue.clear();
    }

    /// `peek()` — the least element, or `None` if empty.
    pub fn peek(&self) -> Option<&T> {
        self.queue.first()
    }

    /// `offer(E e)` / `add(E e)` — insert, then `siftUp`.
    pub fn offer(&mut self, e: T) {
        let i = self.queue.len();
        self.queue.push(e);
        if i != 0 {
            self.sift_up(i);
        }
    }

    /// `add(E e)`.
    pub fn add(&mut self, e: T) {
        self.offer(e);
    }

    /// `poll()` — remove and return the least element, or `None` if empty.
    pub fn poll(&mut self) -> Option<T> {
        let n = self.queue.len();
        if n == 0 {
            return None;
        }
        // swap_remove moves the last element to index 0 (Java's `x = queue[--size]`);
        // siftDown unless that emptied the heap (Java's `if (s != 0)`).
        let result = self.queue.swap_remove(0);
        if !self.queue.is_empty() {
            self.sift_down(0);
        }
        Some(result)
    }

    /// `siftUp(int k, E x)` (comparable form): bubble element at `k` toward the root while it
    /// is strictly less than its parent.
    fn sift_up(&mut self, mut k: usize) {
        while k > 0 {
            let parent = (k - 1) >> 1;
            if self.queue[k].cmp(&self.queue[parent]) != Ordering::Less {
                break;
            }
            self.queue.swap(k, parent);
            k = parent;
        }
    }

    /// `siftDown(int k, E x)` (comparable form): push element at `k` toward the leaves,
    /// preferring the smaller child and the left child on ties (`compareTo > 0` selects right).
    fn sift_down(&mut self, mut k: usize) {
        let n = self.queue.len();
        let half = n >> 1;
        while k < half {
            let mut child = (k << 1) + 1;
            let right = child + 1;
            if right < n && self.queue[child].cmp(&self.queue[right]) == Ordering::Greater {
                child = right;
            }
            if self.queue[k].cmp(&self.queue[child]) != Ordering::Greater {
                break;
            }
            self.queue.swap(k, child);
            k = child;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_sequence_for_seed_0() {
        // Values verified against java.util.Random(0).
        let mut r = Random::new(0);
        assert_eq!(r.next_int(), -1155484576);
        assert_eq!(r.next_int(), -723955400);
        assert_eq!(r.next_int(), 1033096058);
    }

    #[test]
    fn next_int_bound_power_of_two_and_general() {
        let mut r = Random::new(42);
        for _ in 0..1000 {
            let v = r.next_int_bound(100);
            assert!((0..100).contains(&v));
        }
        let mut r2 = Random::new(42);
        for _ in 0..1000 {
            let v = r2.next_int_bound(64);
            assert!((0..64).contains(&v));
        }
    }

    #[test]
    fn next_double_in_unit_interval() {
        let mut r = Random::new(123456789);
        for _ in 0..1000 {
            let d = r.next_double();
            assert!((0.0..1.0).contains(&d));
        }
    }

    #[test]
    fn priority_queue_polls_in_ascending_order() {
        let mut pq = PriorityQueue::new(4);
        for v in [5, 1, 3, 8, 2, 7, 4, 6, 0, 9] {
            pq.offer(v);
        }
        assert_eq!(pq.size(), 10);
        assert_eq!(pq.peek().copied(), Some(0));
        let mut out = Vec::new();
        while let Some(v) = pq.poll() {
            out.push(v);
        }
        assert_eq!(out, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert!(pq.is_empty());
    }

    #[test]
    fn priority_queue_internal_array_matches_jdk() {
        // Heap array layout verified against java.util.PriorityQueue<Integer>.
        // Offering 3,1,4,1,5,9,2,6 yields internal array [1,1,2,3,5,9,4,6].
        let mut pq = PriorityQueue::new(8);
        for v in [3, 1, 4, 1, 5, 9, 2, 6] {
            pq.offer(v);
        }
        assert_eq!(pq.queue, vec![1, 1, 2, 3, 5, 9, 4, 6]);
        // After one poll the JDK array is [1,3,2,6,5,9,4].
        assert_eq!(pq.poll(), Some(1));
        assert_eq!(pq.queue, vec![1, 3, 2, 6, 5, 9, 4]);
    }
}
