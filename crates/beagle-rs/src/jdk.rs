//! Faithful reimplementations of JDK standard-library classes that Beagle uses
//! directly (not part of any Beagle package). Determinism of the port depends on
//! these matching the JDK bit-for-bit.

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
}
