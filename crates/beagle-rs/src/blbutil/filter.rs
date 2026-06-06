//! Port of `blbutil/Filter.java` — accept/reject predicate with include/exclude
//! factories. Java's null checks are not applicable in Rust.

use std::collections::HashSet;
use std::hash::Hash;

/// Port of `blbutil/Filter.java`.
pub enum Filter<E: Eq + Hash> {
    /// `acceptAllFilter()`.
    AcceptAll,
    /// `includeFilter(include)` — accept iff contained.
    Include(HashSet<E>),
    /// `excludeFilter(exclude)` — accept iff not contained.
    Exclude(HashSet<E>),
}

impl<E: Eq + Hash> Filter<E> {
    /// `Filter.acceptAllFilter()`.
    pub fn accept_all() -> Self {
        Filter::AcceptAll
    }

    /// `Filter.includeFilter(Collection)`.
    pub fn include<I: IntoIterator<Item = E>>(include: I) -> Self {
        Filter::Include(include.into_iter().collect())
    }

    /// `Filter.excludeFilter(Collection)`.
    pub fn exclude<I: IntoIterator<Item = E>>(exclude: I) -> Self {
        Filter::Exclude(exclude.into_iter().collect())
    }

    /// `accept(E e)`.
    pub fn accept(&self, e: &E) -> bool {
        match self {
            Filter::AcceptAll => true,
            Filter::Include(s) => s.contains(e),
            Filter::Exclude(s) => !s.contains(e),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filters() {
        let all: Filter<String> = Filter::accept_all();
        assert!(all.accept(&"anything".to_string()));

        let inc = Filter::include(["a".to_string(), "b".to_string()]);
        assert!(inc.accept(&"a".to_string()));
        assert!(!inc.accept(&"c".to_string()));

        let exc = Filter::exclude(["a".to_string(), "b".to_string()]);
        assert!(!exc.accept(&"a".to_string()));
        assert!(exc.accept(&"c".to_string()));
    }
}
