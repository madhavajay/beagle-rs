//! Port of `blbutil/Validate.java` — static helpers for validating command-line arguments.
//!
//! Java's `Map<String,String>` of `key=value` args is a `HashMap<String, String>` here; the
//! `*_arg` helpers `remove` the key (so leftover keys can be detected by `confirm_empty_map`)
//! and validate bounds. Argument errors panic with the Java message (the Java
//! `IllegalArgumentException`); `confirm_empty_map`/`get_file` terminate via `Utilities::exit`
//! as in Java.

use std::collections::HashMap;
use std::path::PathBuf;

use super::{consts, Utilities};

/// `argsToMap(String[] args, char delim)`.
pub fn args_to_map(args: &[String], delim: char) -> HashMap<String, String> {
    let mut arg_map = HashMap::new();
    for arg in args {
        match arg.find(delim) {
            Some(index) => {
                assert!(index != 0, "missing key in key-value pair: {arg}");
                assert!(
                    index != arg.len() - delim.len_utf8(),
                    "missing value in key-value pair: {arg}"
                );
                let key = arg[..index].to_string();
                let value = arg[index + delim.len_utf8()..].to_string();
                assert!(!arg_map.contains_key(&key), "duplicate arguments: {key}");
                arg_map.insert(key, value);
            }
            None => panic!("missing delimiter character ({delim}): {arg}"),
        }
    }
    arg_map
}

/// `confirmEmptyMap(Map<String,String> argsMap)` — exits if any keys remain.
pub fn confirm_empty_map(args_map: &HashMap<String, String>) {
    if !args_map.is_empty() {
        let mut sb = String::from("Error: unrecognized parameter");
        sb.push_str(if args_map.len() == 1 { ":" } else { "s:" });
        // Java iterates a HashMap keySet (unspecified order); message order is not
        // output-significant. Sort for determinism.
        let mut keys: Vec<&String> = args_map.keys().collect();
        keys.sort();
        for key in keys {
            sb.push(' ');
            sb.push_str(key);
            sb.push('=');
            sb.push_str(&args_map[key]);
        }
        Utilities::exit(&sb);
    }
}

/// `getFile(String filename)` — `None` if `filename` is `None`; exits if the file is empty,
/// missing, or a directory.
pub fn get_file(filename: Option<&str>) -> Option<PathBuf> {
    let filename = filename?;
    assert!(!filename.is_empty(), "filename is empty string");
    let path = PathBuf::from(filename);
    let err = if !path.exists() {
        Some("File does not exist")
    } else if path.is_dir() {
        Some("File cannot be a directory")
    } else {
        None
    };
    if let Some(err) = err {
        Utilities::exit(&format!(
            "{nl}Error     :  {err}{nl}Filename  :  {}",
            path.display(),
            nl = consts::NL
        ));
    }
    Some(path)
}

/// `intArg(...)`.
pub fn int_arg(
    key: &str,
    map: &mut HashMap<String, String>,
    is_required: bool,
    default_value: i32,
    min: i32,
    max: i32,
) -> i32 {
    check_int_value(key, default_value, min, max);
    match map.remove(key) {
        None => {
            assert!(!is_required, "missing {key} argument");
            default_value
        }
        Some(value) => parse_int(key, &value, min, max),
    }
}

/// `longArg(...)`.
pub fn long_arg(
    key: &str,
    map: &mut HashMap<String, String>,
    is_required: bool,
    default_value: i64,
    min: i64,
    max: i64,
) -> i64 {
    check_long_value(key, default_value, min, max);
    match map.remove(key) {
        None => {
            assert!(!is_required, "missing {key} argument");
            default_value
        }
        Some(value) => parse_long(key, &value, min, max),
    }
}

/// `floatArg(...)`.
pub fn float_arg(
    key: &str,
    map: &mut HashMap<String, String>,
    is_required: bool,
    default_value: f32,
    min: f32,
    max: f32,
) -> f32 {
    check_float_value(key, default_value, min, max);
    match map.remove(key) {
        None => {
            assert!(!is_required, "missing {key} argument");
            default_value
        }
        Some(value) => parse_float(key, &value, min, max),
    }
}

/// `doubleArg(...)`.
pub fn double_arg(
    key: &str,
    map: &mut HashMap<String, String>,
    is_required: bool,
    default_value: f64,
    min: f64,
    max: f64,
) -> f64 {
    check_double_value(key, default_value, min, max);
    match map.remove(key) {
        None => {
            assert!(!is_required, "missing {key} argument");
            default_value
        }
        Some(value) => parse_double(key, &value, min, max),
    }
}

/// `booleanArg(...)` — accepts `true/t/false/f` (case-insensitive).
pub fn boolean_arg(
    key: &str,
    map: &mut HashMap<String, String>,
    is_required: bool,
    default_value: bool,
) -> bool {
    match map.remove(key) {
        None => {
            assert!(!is_required, "missing {key} argument");
            default_value
        }
        Some(value) => parse_boolean(&value),
    }
}

/// `stringArg(...)` — value may be `None`. `possible_values = None` means any non-empty
/// string (and `None`) is allowed.
pub fn string_arg(
    key: &str,
    map: &mut HashMap<String, String>,
    is_required: bool,
    default_value: Option<&str>,
    possible_values: Option<&[&str]>,
) -> Option<String> {
    check_string_value(key, default_value, possible_values);
    match map.remove(key) {
        None => {
            assert!(!is_required, "missing {key} argument");
            default_value.map(str::to_string)
        }
        Some(value) => {
            check_string_value(key, Some(&value), possible_values);
            Some(value)
        }
    }
}

fn parse_int(key: &str, to_parse: &str, min: i32, max: i32) -> i32 {
    let i: i32 = to_parse
        .parse()
        .unwrap_or_else(|_| panic!("{to_parse} is not a number"));
    check_int_value(key, i, min, max);
    i
}

fn parse_long(key: &str, to_parse: &str, min: i64, max: i64) -> i64 {
    let l: i64 = to_parse
        .parse()
        .unwrap_or_else(|_| panic!("{to_parse} is not a number"));
    check_long_value(key, l, min, max);
    l
}

fn parse_float(key: &str, to_parse: &str, min: f32, max: f32) -> f32 {
    let f: f32 = to_parse
        .parse()
        .unwrap_or_else(|_| panic!("{to_parse} is not a number"));
    check_float_value(key, f, min, max);
    f
}

fn parse_double(key: &str, to_parse: &str, min: f64, max: f64) -> f64 {
    let d: f64 = to_parse
        .parse()
        .unwrap_or_else(|_| panic!("{to_parse} is not a number"));
    check_double_value(key, d, min, max);
    d
}

fn parse_boolean(s: &str) -> bool {
    if s.eq_ignore_ascii_case("true") || s.eq_ignore_ascii_case("t") {
        true
    } else if s.eq_ignore_ascii_case("false") || s.eq_ignore_ascii_case("f") {
        false
    } else {
        panic!("{s} is not \"true\" or \"false\"");
    }
}

fn check_int_value(key: &str, value: i32, min: i32, max: i32) {
    let s = if min > max {
        Some(format!("min={min} > max={max}"))
    } else if value < min {
        Some(format!("value={value} < {min}"))
    } else if value > max {
        Some(format!("value={value} > {max}"))
    } else {
        None
    };
    if let Some(s) = s {
        panic!("Error in \"{key}\" argument: {s}");
    }
}

fn check_long_value(key: &str, value: i64, min: i64, max: i64) {
    let s = if min > max {
        Some(format!("min={min} > max={max}"))
    } else if value < min {
        Some(format!("value={value} < {min}"))
    } else if value > max {
        Some(format!("value={value} > {max}"))
    } else {
        None
    };
    if let Some(s) = s {
        panic!("Error in \"{key}\" argument: {s}");
    }
}

fn check_float_value(key: &str, value: f32, min: f32, max: f32) {
    let s = if value.is_nan() {
        Some(format!("value={value}"))
    } else if min > max {
        Some(format!("min={min} > max={max}"))
    } else if value < min {
        Some(format!("value={value} < {min}"))
    } else if value > max {
        Some(format!("value={value} > {max}"))
    } else {
        None
    };
    if let Some(s) = s {
        panic!("Error in \"{key}\" argument: {s}");
    }
}

fn check_double_value(key: &str, value: f64, min: f64, max: f64) {
    let s = if value.is_nan() {
        Some(format!("value={value}"))
    } else if min > max {
        Some(format!("min={min} > max={max}"))
    } else if value < min {
        Some(format!("value={value} < {min}"))
    } else if value > max {
        Some(format!("value={value} > {max}"))
    } else {
        None
    };
    if let Some(s) = s {
        panic!("Error in \"{key}\" argument: {s}");
    }
}

fn check_string_value(key: &str, value: Option<&str>, possible_values: Option<&[&str]>) {
    if let Some(possible) = possible_values {
        let found = possible
            .iter()
            .any(|s| value.is_some_and(|v| s.eq_ignore_ascii_case(v)));
        if !found {
            panic!(
                "Error in \"{key}\" argument: \"{}\" is not in {possible:?}",
                value.unwrap_or("null")
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn map(pairs: &[(&str, &str)]) -> HashMap<String, String> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    #[test]
    fn args_to_map_parses_pairs() {
        let m = args_to_map(&["gt=in.vcf".into(), "nthreads=4".into()], '=');
        assert_eq!(m["gt"], "in.vcf");
        assert_eq!(m["nthreads"], "4");
    }

    #[test]
    fn int_arg_default_and_value_and_remove() {
        let mut m = map(&[("window", "40")]);
        assert_eq!(int_arg("window", &mut m, false, 100, 1, 1000), 40);
        assert!(!m.contains_key("window")); // removed
                                            // default when absent
        assert_eq!(int_arg("overlap", &mut m, false, 2, 0, 10), 2);
    }

    #[test]
    #[should_panic(expected = "Error in \"window\" argument: value=2000 > 1000")]
    fn int_arg_out_of_range() {
        let mut m = map(&[("window", "2000")]);
        int_arg("window", &mut m, false, 100, 1, 1000);
    }

    #[test]
    fn boolean_and_string_and_float() {
        let mut m = map(&[("gp", "TRUE"), ("err", "1.0e-4"), ("out", "result")]);
        assert!(boolean_arg("gp", &mut m, false, false));
        assert!((float_arg("err", &mut m, false, 1e-4, 0.0, 1.0) - 1e-4).abs() < 1e-12);
        assert_eq!(
            string_arg("out", &mut m, true, None, None),
            Some("result".to_string())
        );
        assert!(m.is_empty());
    }

    #[test]
    fn string_arg_possible_values() {
        let mut m = map(&[("mode", "B")]);
        assert_eq!(
            string_arg("mode", &mut m, false, Some("a"), Some(&["a", "b", "c"])),
            Some("B".to_string())
        );
    }

    #[test]
    fn confirm_empty_map_ok_when_empty() {
        confirm_empty_map(&HashMap::new());
    }
}
