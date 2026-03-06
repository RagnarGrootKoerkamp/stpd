use std::{cmp::Ordering, ops::Range};

/// Implementation of lex-min STPD.
pub struct Stpd {
    /// The input text.
    // TODO: Relative LZ encoding.
    text: Vec<u8>,
    /// STPD samples, covering all left-most occurrences of all right-maximal extensions.
    /// Sorted in co-lex order.
    /// spa: sampled prefix array
    // TODO: BTree instead so we can insert samples.
    spa: Vec<usize>,
}
impl Stpd {
    pub fn new(text: &[u8]) -> Self {
        todo!();
    }

    /// Find the range of `spa` that has `q` as a suffix.
    // TODO: Do everything in reverse instead, so we have normal lex comparisons?
    fn binary_search(&self, q: &[u8]) -> Range<usize> {
        let start = self
            .spa
            .binary_search_by(|idx| match cmp_colex(&self.text[..*idx], q) {
                Ordering::Equal => Ordering::Greater,
                x => x,
            })
            .unwrap_err();
        let end = self
            .spa
            .binary_search_by(|idx| match cmp_colex(&self.text[..*idx], q) {
                Ordering::Equal => Ordering::Less,
                x => x,
            })
            .unwrap_err();
        start..end
    }

    /// Return the leftmost occurrence of `q` in `text`, if it exists.
    fn search_rme(&self, q: &[u8]) -> Option<usize> {
        assert!(!q.is_empty());
        let range = self.binary_search(q);
        if range.is_empty() {
            return None;
        }
        // Find the smallest index in the range.
        let rme_idx = *self.spa[range].iter().min().unwrap();
        Some(rme_idx - (q.len() - 1))
    }

    pub fn locate_one(&self, q: &[u8]) -> Option<usize> {
        // Number of matched characters of q.
        let mut i = 0;
        // Start position of leftmost occurrence in text of q[..i].
        let mut pos = 0;
        while pos < q.len() {
            if self.text[pos + i] == q[i] {
                i += 1;
                continue;
            }
            // q[..=i] does not occur at `pos`, so is either an RME or does not occur at all.
            pos = self.search_rme(&q[..=i])?;
        }
        Some(pos)
    }
}

/// Length of longest common suffix.
fn lcs(a: &[u8], b: &[u8]) -> usize {
    let mut i = 0;
    while i < Ord::min(a.len(), b.len()) && a[a.len() - 1 - i] == b[b.len() - 1 - i] {
        i += 1;
    }
    return i;
}

/// co-lex compare q with a text prefix.
/// Returns `Equal` when `q` is a suffix of `text`.
fn cmp_colex(text: &[u8], q: &[u8]) -> Ordering {
    let l = lcs(text, q);
    if l == q.len() {
        return Ordering::Equal;
    }
    if l == text.len() {
        return Ordering::Less;
    }
    return Ord::cmp(&text[text.len() - 1 - l], &q[q.len() - 1 - l]);
}
