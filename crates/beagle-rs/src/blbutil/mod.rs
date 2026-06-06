//! Port of the Java `blbutil` package: base utility types shared across Beagle.
//!
//! Ported so far: `Const` (constants), `FloatArray`, `DoubleArray`, `FloatList`,
//! `BitArray`, `StringUtil`. Remaining (IO/BGZIP/validation) land in later chunks.

mod bgzip;
mod bit_array;
mod filter;
mod float_arrays;
mod string_util;
mod utilities;

pub use bgzip::{write_empty_block, BgzipOutputStream, MAX_INPUT_BYTES};
pub use bit_array::BitArray;
pub use filter::Filter;
pub use float_arrays::{DoubleArray, FloatArray, FloatList};
pub use string_util::StringUtil;
pub use utilities::Utilities;

/// Port of `blbutil/Const.java` — string/character/numeric constants. (Java field
/// names are lower-case, e.g. `Const.tab`; here they follow Rust `SCREAMING_CASE`.)
pub mod consts {
    /// `Const.nl` — the platform line separator. The byte-for-byte reference runs on
    /// Linux, where `System.getProperty("line.separator")` is `"\n"`.
    pub const NL: &str = "\n";
    /// `Const.MISSING_DATA_STRING`.
    pub const MISSING_DATA_STRING: &str = ".";
    /// `Const.MISSING_DATA_CHAR`.
    pub const MISSING_DATA_CHAR: char = '.';
    /// `Const.colon`.
    pub const COLON: char = ':';
    /// `Const.hyphen`.
    pub const HYPHEN: char = '-';
    /// `Const.tab`.
    pub const TAB: char = '\t';
    /// `Const.semicolon`.
    pub const SEMICOLON: char = ';';
    /// `Const.comma`.
    pub const COMMA: char = ',';
    /// `Const.phasedSep`.
    pub const PHASED_SEP: char = '|';
    /// `Const.unphasedSep`.
    pub const UNPHASED_SEP: char = '/';
    /// `Const.giga` = 1,000,000,000.
    pub const GIGA: i32 = 1_000_000_000;
    /// `Const.mega` = 1,000,000.
    pub const MEGA: i32 = 1_000_000;
}

// ---- Java `long` shift semantics (the JVM masks the shift count to 6 bits) ----

/// Java `x << n` on `long`.
#[inline]
pub(crate) fn jshl_l(x: i64, n: i32) -> i64 {
    x.wrapping_shl(n as u32)
}

/// Java `x >> n` on `long` (arithmetic).
#[inline]
pub(crate) fn jshr_l(x: i64, n: i32) -> i64 {
    x.wrapping_shr(n as u32)
}

/// Java `x >>> n` on `long` (logical).
#[inline]
pub(crate) fn jushr_l(x: i64, n: i32) -> i64 {
    (x as u64).wrapping_shr(n as u32) as i64
}

/// `Float.floatToIntBits` (collapses all NaN to the canonical `0x7fc00000`).
#[inline]
pub(crate) fn float_to_int_bits(f: f32) -> i32 {
    if f.is_nan() {
        0x7fc0_0000u32 as i32
    } else {
        f.to_bits() as i32
    }
}

/// `Double.doubleToLongBits` (collapses all NaN to the canonical `0x7ff8…`).
#[inline]
pub(crate) fn double_to_long_bits(d: f64) -> i64 {
    if d.is_nan() {
        0x7ff8_0000_0000_0000u64 as i64
    } else {
        d.to_bits() as i64
    }
}

/// `java.util.Arrays.binarySearch(float[], fromIndex, toIndex, key)` semantics
/// (NaN treated as greatest; `-0.0 < 0.0`).
pub(crate) fn binary_search_f32(a: &[f32], from: i32, to: i32, key: f32) -> i32 {
    let mut low = from;
    let mut high = to - 1;
    while low <= high {
        let mid = ((low as u32).wrapping_add(high as u32) >> 1) as i32;
        let mid_val = a[mid as usize];
        if mid_val < key {
            low = mid + 1;
        } else if mid_val > key {
            high = mid - 1;
        } else {
            let mid_bits = float_to_int_bits(mid_val);
            let key_bits = float_to_int_bits(key);
            if mid_bits == key_bits {
                return mid;
            } else if mid_bits < key_bits {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
    }
    -(low + 1)
}

/// `java.util.Arrays.binarySearch(double[], fromIndex, toIndex, key)` semantics.
pub(crate) fn binary_search_f64(a: &[f64], from: i32, to: i32, key: f64) -> i32 {
    let mut low = from;
    let mut high = to - 1;
    while low <= high {
        let mid = ((low as u32).wrapping_add(high as u32) >> 1) as i32;
        let mid_val = a[mid as usize];
        if mid_val < key {
            low = mid + 1;
        } else if mid_val > key {
            high = mid - 1;
        } else {
            let mid_bits = double_to_long_bits(mid_val);
            let key_bits = double_to_long_bits(key);
            if mid_bits == key_bits {
                return mid;
            } else if mid_bits < key_bits {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
    }
    -(low + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn long_shifts_match_java() {
        assert_eq!(jshl_l(1, 63), i64::MIN);
        assert_eq!(jshl_l(1, 64), 1); // shift count masked to 6 bits
        assert_eq!(jshr_l(i64::MIN, 63), -1); // arithmetic
        assert_eq!(jushr_l(i64::MIN, 63), 1); // logical
        assert_eq!(jushr_l(-1, 32), 0xffff_ffff);
    }

    #[test]
    fn float_bits_canonical_nan() {
        assert_eq!(float_to_int_bits(f32::NAN), 0x7fc0_0000u32 as i32);
        assert_eq!(float_to_int_bits(0.0), 0);
        assert_eq!(float_to_int_bits(-0.0), i32::MIN);
        assert_eq!(
            double_to_long_bits(f64::NAN),
            0x7ff8_0000_0000_0000u64 as i64
        );
    }

    #[test]
    fn binary_search_floats() {
        let a = [1.0f32, 3.0, 5.0, 7.0];
        assert_eq!(binary_search_f32(&a, 0, 4, 5.0), 2);
        assert_eq!(binary_search_f32(&a, 0, 4, 4.0), -3);
        assert_eq!(binary_search_f32(&a, 0, 4, 0.0), -1);
        let d = [1.0f64, 2.5, 9.0];
        assert_eq!(binary_search_f64(&d, 0, 3, 2.5), 1);
        assert_eq!(binary_search_f64(&d, 0, 3, 3.0), -3);
    }
}
