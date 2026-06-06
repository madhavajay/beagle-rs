//! Cross-language parity for `java.util.Random` + `Utilities.shuffle`. The Rust port
//! must reproduce `fixtures/random/parity.txt` (from `tools/java/RandomParityDriver.java`)
//! exactly. This is the determinism foundation for byte-for-byte phasing/imputation.

use beagle_rs::blbutil::Utilities;
use beagle_rs::jdk::Random;

fn csv<T: ToString>(v: &[T]) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn arrays_to_string(a: &[i32]) -> String {
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
    out.push_str("# random parity (java.util.Random)\n");

    let mut r = Random::new(42);
    let ni: Vec<i32> = (0..8).map(|_| r.next_int()).collect();
    out.push_str(&format!("nextInt={}\n", csv(&ni)));
    let nb100: Vec<i32> = (0..8).map(|_| r.next_int_bound(100)).collect();
    out.push_str(&format!("nextIntBound100={}\n", csv(&nb100)));
    let nb64: Vec<i32> = (0..8).map(|_| r.next_int_bound(64)).collect();
    out.push_str(&format!("nextIntBound64={}\n", csv(&nb64)));
    let nl: Vec<i64> = (0..4).map(|_| r.next_long()).collect();
    out.push_str(&format!("nextLong={}\n", csv(&nl)));
    let nd: Vec<i64> = (0..4).map(|_| r.next_double().to_bits() as i64).collect();
    out.push_str(&format!("nextDoubleBits={}\n", csv(&nd)));
    let nbool: Vec<i32> = (0..8).map(|_| i32::from(r.next_boolean())).collect();
    out.push_str(&format!("nextBoolean={}\n", csv(&nbool)));

    let mut ia: Vec<i32> = (0..20).collect();
    let mut rs = Random::new(7);
    Utilities::shuffle(&mut ia, &mut rs);
    out.push_str(&format!("shuffle={}\n", arrays_to_string(&ia)));

    out
}

#[test]
fn random_parity_matches_java_reference() {
    let expected = include_str!("../../../fixtures/random/parity.txt");
    assert_eq!(
        transcript(),
        expected,
        "Rust java.util.Random port diverged from the JDK"
    );
}
