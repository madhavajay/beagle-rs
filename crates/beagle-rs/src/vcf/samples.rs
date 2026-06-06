//! Port of `vcf/Samples.java` — an immutable ordered list of samples, each with an id
//! and a haploid/diploid flag.

/// Port of `vcf/Samples.java`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Samples {
    ids: Vec<String>,
    is_diploid: Vec<bool>,
}

fn check_for_nulls_and_duplicates(ids: &[String]) {
    let mut sorted: Vec<&String> = ids.iter().collect();
    sorted.sort();
    if let Some(first) = sorted.first() {
        if first.is_empty() {
            panic!("Empty string identifier");
        }
    }
    for j in 1..sorted.len() {
        if sorted[j].is_empty() {
            panic!("Empty string identifier");
        }
        if sorted[j] == sorted[j - 1] {
            eprintln!("Warning: duplicate sample identifier: {}", sorted[j]);
        }
    }
}

impl Samples {
    /// `new Samples(String[] ids, boolean[] isDiploid)`.
    pub fn new(ids: &[String], is_diploid: &[bool]) -> Self {
        assert!(
            ids.len() == is_diploid.len(),
            "ids.length ({}) != isDiploid.length ({})",
            ids.len(),
            is_diploid.len()
        );
        check_for_nulls_and_duplicates(ids);
        Samples {
            ids: ids.to_vec(),
            is_diploid: is_diploid.to_vec(),
        }
    }

    /// `Samples.combine(first, second)` — concatenates the two lists in order.
    pub fn combine(first: &Samples, second: &Samples) -> Samples {
        let mut ids = Vec::with_capacity(first.ids.len() + second.ids.len());
        ids.extend_from_slice(&first.ids);
        ids.extend_from_slice(&second.ids);
        let mut is_diploid = Vec::with_capacity(ids.len());
        is_diploid.extend_from_slice(&first.is_diploid);
        is_diploid.extend_from_slice(&second.is_diploid);
        Samples::new(&ids, &is_diploid)
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.ids.len() as i32
    }

    /// `id(int index)`.
    pub fn id(&self, index: i32) -> &str {
        &self.ids[index as usize]
    }

    /// `ids()`.
    pub fn ids(&self) -> Vec<String> {
        self.ids.clone()
    }

    /// `isDiploid(int sample)`.
    pub fn is_diploid(&self, sample: i32) -> bool {
        self.is_diploid[sample as usize]
    }
}

impl std::fmt::Display for Samples {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // java.util.Arrays.toString(String[]): "[A, B, C]"
        f.write_str("[")?;
        for (i, id) in self.ids.iter().enumerate() {
            if i > 0 {
                f.write_str(", ")?;
            }
            f.write_str(id)?;
        }
        f.write_str("]")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn basics() {
        let samp = Samples::new(&s(&["A", "B", "C"]), &[true, false, true]);
        assert_eq!(samp.size(), 3);
        assert_eq!(samp.id(1), "B");
        assert!(samp.is_diploid(0));
        assert!(!samp.is_diploid(1));
        assert_eq!(samp.ids(), s(&["A", "B", "C"]));
        assert_eq!(samp.to_string(), "[A, B, C]");
    }

    #[test]
    fn combine_concatenates() {
        let a = Samples::new(&s(&["A", "B"]), &[true, true]);
        let b = Samples::new(&s(&["C"]), &[false]);
        let c = Samples::combine(&a, &b);
        assert_eq!(c.ids(), s(&["A", "B", "C"]));
        assert!(!c.is_diploid(2));
        assert_eq!(a, Samples::new(&s(&["A", "B"]), &[true, true]));
    }

    #[test]
    #[should_panic]
    fn rejects_length_mismatch() {
        let _ = Samples::new(&s(&["A", "B"]), &[true]);
    }

    #[test]
    #[should_panic]
    fn rejects_empty_id() {
        let _ = Samples::new(&s(&["A", ""]), &[true, true]);
    }
}
