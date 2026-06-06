//! Port of `blbutil/StringUtil.java` — field counting/splitting on a delimiter or on
//! white space.
//!
//! Java operates on UTF-16 `char` units with `c <= ' '` for white space. This port
//! operates on bytes, which is byte-identical for the ASCII delimiters Beagle uses
//! (tab, etc.) and for white space (all chars `<= U+0020` are single-byte ASCII).
//! Splitting on an ASCII delimiter never lands inside a multi-byte UTF-8 sequence, so
//! the returned `&str` slices are always valid.

/// Port of `blbutil/StringUtil.java` (static utility methods).
pub struct StringUtil;

fn trim_ws(s: &str) -> &str {
    let b = s.as_bytes();
    let mut start = 0usize;
    let mut end = b.len();
    while start < end && b[start] <= b' ' {
        start += 1;
    }
    while end > start && b[end - 1] <= b' ' {
        end -= 1;
    }
    &s[start..end]
}

impl StringUtil {
    /// `countFields(String s, char delimiter)` — number of delimited fields (0-length
    /// string still yields 1, matching Java).
    pub fn count_fields(s: &str, delimiter: u8) -> i32 {
        let cnt = s.as_bytes().iter().filter(|&&b| b == delimiter).count() as i32;
        cnt + 1
    }

    /// `countFields(String s, char delimiter, int max)`.
    pub fn count_fields_max(s: &str, delimiter: u8, max: i32) -> i32 {
        let mut cnt = 0i32;
        let max_cnt = max - 1;
        for &b in s.as_bytes() {
            if cnt >= max_cnt {
                break;
            }
            if b == delimiter {
                cnt += 1;
            }
        }
        (cnt + 1).min(max)
    }

    /// `getFields(String s, char delimiter)`.
    pub fn get_fields(s: &str, delimiter: u8) -> Vec<&str> {
        let bytes = s.as_bytes();
        let n = Self::count_fields(s, delimiter);
        let mut fields = Vec::with_capacity(n as usize);
        let mut start = 0usize;
        for _ in 0..n {
            match bytes[start..].iter().position(|&b| b == delimiter) {
                Some(off) => {
                    let end = start + off;
                    fields.push(&s[start..end]);
                    start = end + 1;
                }
                None => fields.push(&s[start..]),
            }
        }
        fields
    }

    /// `getFields(String s, char delimiter, int limit)`.
    pub fn get_fields_limit(s: &str, delimiter: u8, limit: i32) -> Vec<&str> {
        assert!(limit >= 2, "limit: {}", limit);
        let count = Self::count_fields_max(s, delimiter, limit);
        let mut fields = Vec::with_capacity(count.max(0) as usize);
        if count > 0 {
            let bytes = s.as_bytes();
            let mut start = 0usize;
            for _ in 0..(count - 1) {
                let off = bytes[start..]
                    .iter()
                    .position(|&b| b == delimiter)
                    .expect("delimiter expected within counted fields");
                let end = start + off;
                fields.push(&s[start..end]);
                start = end + 1;
            }
            fields.push(&s[start..]);
        }
        fields
    }

    /// `countFields(String s)` — white-space delimited field count.
    pub fn count_fields_ws(s: &str) -> i32 {
        let mut cnt = 0i32;
        let mut in_field = false;
        for &b in s.as_bytes() {
            if b > b' ' {
                if !in_field {
                    cnt += 1;
                    in_field = true;
                }
            } else {
                in_field = false;
            }
        }
        cnt
    }

    /// `getFields(String s)` — trim and split on white space.
    pub fn get_fields_ws(s: &str) -> Vec<&str> {
        let bytes = s.as_bytes();
        let mut out = Vec::new();
        let mut i = 0usize;
        while i < bytes.len() {
            if bytes[i] > b' ' {
                let start = i;
                while i < bytes.len() && bytes[i] > b' ' {
                    i += 1;
                }
                out.push(&s[start..i]);
            } else {
                i += 1;
            }
        }
        out
    }

    /// `getFields(String s, int limit)` — trim, then split on the first `limit-1`
    /// white-space delimiters (the last field is the remainder).
    pub fn get_fields_ws_limit(s: &str, limit: i32) -> Vec<&str> {
        assert!(limit >= 2, "limit: {}", limit);
        let t = trim_ws(s);
        let bytes = t.as_bytes();
        let mut out = Vec::new();
        let mut i = 0usize;
        let mut idx = 0i32;
        while i < bytes.len() {
            while i < bytes.len() && bytes[i] <= b' ' {
                i += 1;
            }
            if i >= bytes.len() {
                break;
            }
            let start = i;
            while i < bytes.len() && bytes[i] > b' ' {
                i += 1;
            }
            if idx < limit - 1 {
                out.push(&t[start..i]);
                idx += 1;
            } else {
                out.push(&t[start..]);
                break;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAB: u8 = b'\t';

    #[test]
    fn count_and_get_fields_delim() {
        assert_eq!(StringUtil::count_fields("a\tb\tc", TAB), 3);
        assert_eq!(StringUtil::count_fields("", TAB), 1);
        assert_eq!(StringUtil::count_fields("a\t", TAB), 2);
        assert_eq!(StringUtil::get_fields("a\tb\tc", TAB), vec!["a", "b", "c"]);
        assert_eq!(StringUtil::get_fields("a\t\tc", TAB), vec!["a", "", "c"]);
        assert_eq!(StringUtil::get_fields("a\t", TAB), vec!["a", ""]);
        assert_eq!(StringUtil::get_fields("", TAB), vec![""]);
    }

    #[test]
    fn count_and_get_fields_limit() {
        assert_eq!(StringUtil::count_fields_max("a\tb\tc\td", TAB, 2), 2);
        assert_eq!(
            StringUtil::get_fields_limit("a\tb\tc\td", TAB, 2),
            vec!["a", "b\tc\td"]
        );
        assert_eq!(
            StringUtil::get_fields_limit("a\tb\tc", TAB, 10),
            vec!["a", "b", "c"]
        );
    }

    #[test]
    fn whitespace_fields() {
        assert_eq!(StringUtil::count_fields_ws("  a  bb   ccc  "), 3);
        assert_eq!(
            StringUtil::get_fields_ws("  a  bb   ccc  "),
            vec!["a", "bb", "ccc"]
        );
        assert_eq!(StringUtil::get_fields_ws("   "), Vec::<&str>::new());
        assert_eq!(StringUtil::get_fields_ws("solo"), vec!["solo"]);
    }

    #[test]
    fn whitespace_fields_limit() {
        assert_eq!(
            StringUtil::get_fields_ws_limit("  a  b  c  d  ", 2),
            vec!["a", "b  c  d"]
        );
        assert_eq!(
            StringUtil::get_fields_ws_limit("a b c", 10),
            vec!["a", "b", "c"]
        );
        assert_eq!(
            StringUtil::get_fields_ws_limit("   ", 3),
            Vec::<&str>::new()
        );
    }

    #[test]
    fn utf8_safe_split_on_ascii_delim() {
        // multi-byte chars survive splitting on an ASCII delimiter
        let fields = StringUtil::get_fields("café\tnaïve", TAB);
        assert_eq!(fields, vec!["café", "naïve"]);
    }
}
