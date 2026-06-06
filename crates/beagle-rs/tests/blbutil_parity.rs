//! Cross-language parity for the `blbutil` package: the Rust port must reproduce the
//! transcript emitted by `tools/java/BlbutilParityDriver.java`
//! (`fixtures/blbutil/parity.txt`). Verifies BitArray's byte-exact `long`-word layout
//! and StringUtil's splitting semantics against the Java reference.

use beagle_rs::blbutil::{BitArray, StringUtil};

fn longs_to_string(a: &[i64]) -> String {
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

fn transcript() -> String {
    let mut out = String::new();
    out.push_str("# blbutil parity transcript (Beagle 5.5 27Feb25.75f)\n");

    let mut b = BitArray::new(200);
    for j in 0..200 {
        if j % 7 == 0 || j % 11 == 3 {
            b.set(j);
        }
    }
    out.push_str(&format!("bitarray words={}\n", longs_to_string(&b.to_long_array())));
    out.push_str(&format!(
        "bitarray getAsInt[0,7,63,64,193]={},{},{},{},{}\n",
        b.get_as_int(0),
        b.get_as_int(7),
        b.get_as_int(63),
        b.get_as_int(64),
        b.get_as_int(193)
    ));

    let r = b.restrict(5, 130);
    out.push_str(&format!(
        "bitarray restrict5_130 size={} words={}\n",
        r.size(),
        longs_to_string(&r.to_long_array())
    ));

    out.push_str(&format!(
        "bitarray hash 0_200={} hash 5_130={}\n",
        b.hash(0, 200),
        b.hash(5, 130)
    ));

    let mut d = BitArray::new(200);
    d.copy_from(&b, 10, 150);
    out.push_str(&format!(
        "bitarray copyfrom10_150 words={}\n",
        longs_to_string(&d.to_long_array())
    ));

    out.push_str(&format!(
        "bitarray longHashCode={}\n",
        BitArray::long_hash_code(0x123456789abcdefi64)
    ));

    out.push_str(&format!(
        "split_tab={}\n",
        StringUtil::get_fields("a\tb\t\tc\t", b'\t').join("|")
    ));
    out.push_str(&format!(
        "split_tab_limit2={}\n",
        StringUtil::get_fields_limit("a\tb\tc\td", b'\t', 2).join("|")
    ));
    out.push_str(&format!(
        "count_tab={} count_tab_max2={}\n",
        StringUtil::count_fields("a\tb\tc", b'\t'),
        StringUtil::count_fields_max("a\tb\tc\td", b'\t', 2)
    ));
    out.push_str(&format!(
        "split_ws={}\n",
        StringUtil::get_fields_ws("  hello   world  foo ").join("|")
    ));
    out.push_str(&format!(
        "split_ws_limit2={}\n",
        StringUtil::get_fields_ws_limit("  hello   world  foo ", 2).join("|")
    ));
    out.push_str(&format!(
        "count_ws={}\n",
        StringUtil::count_fields_ws("  hello   world  foo ")
    ));

    out
}

#[test]
fn blbutil_parity_matches_java_reference() {
    let expected = include_str!("../../../fixtures/blbutil/parity.txt");
    assert_eq!(
        transcript(),
        expected,
        "Rust blbutil transcript diverged from the Java reference"
    );
}
