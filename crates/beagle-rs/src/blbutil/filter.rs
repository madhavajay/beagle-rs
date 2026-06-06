//! Port of `blbutil/Filter.java` — accept/reject predicate with include/exclude
//! factories. Java's null checks are not applicable in Rust.

use std::collections::HashSet;
use std::hash::Hash;

/// Port of `blbutil/Filter.java`.
///
/// Java's `Filter<E>` is a `@FunctionalInterface` (`boolean accept(E)`), so callers can
/// supply arbitrary lambdas (e.g. `FilterUtil`'s marker/chrom-interval filters). The common
/// set-membership cases are modeled as enum variants; the `Predicate` variant carries an
/// arbitrary closure for the rest.
pub enum Filter<E: Eq + Hash> {
    /// `acceptAllFilter()`.
    AcceptAll,
    /// `includeFilter(include)` — accept iff contained.
    Include(HashSet<E>),
    /// `excludeFilter(exclude)` — accept iff not contained.
    Exclude(HashSet<E>),
    /// An arbitrary predicate (a Java lambda `Filter<E>`).
    Predicate(Box<dyn Fn(&E) -> bool>),
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

    /// A `Filter` from an arbitrary predicate (a Java lambda `Filter<E>`).
    pub fn predicate<F: Fn(&E) -> bool + 'static>(f: F) -> Self {
        Filter::Predicate(Box::new(f))
    }

    /// `accept(E e)`.
    pub fn accept(&self, e: &E) -> bool {
        match self {
            Filter::AcceptAll => true,
            Filter::Include(s) => s.contains(e),
            Filter::Exclude(s) => !s.contains(e),
            Filter::Predicate(f) => f(e),
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
