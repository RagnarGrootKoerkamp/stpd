use std::{cmp::Ordering, ops::Range};

use itertools::Itertools;

/// Implementation of lex-min STPD.
pub struct Stpd {
    /// The input text.
    // TODO: Relative LZ encoding.
    text: Vec<u8>,
    /// STPD samples, covering all left-most occurrences of all right-maximal extensions.
    /// Sorted in co-lex order.
    /// spa: sampled prefix array
    /// `spa[0]` is the 'trivial' empty anchor at the root of the text.
    // TODO: BTree instead so we can insert samples.
    // TODO: Efficient range min.
    spa: Vec<Anchor>,
}

/// A right maximal extension is a substring sA of T such that s occurs at least twice in the text.
/// We do path compression, which means that the A only sampled if `sA` is
/// leftmost, but `s` has another occurrence further to the left.
/// `pos` indicates the length of the leftmost text prefix `text[..pos]` ending in `sA`.
/// This position is an RME for an interval of length of `s`:
/// - For very short `s'`, `s'A` will already occur to the left of `pos`.
/// - For very long `s''`, this will already be the leftmost occurrence of `s''`.
struct Anchor {
    /// Position in the text of the RME sample.
    /// The prefix `text[..pos]` ends in `sA`.
    pos: usize,
    /// The length of the shortest `sA` such that `sA` is an RME.
    /// This is the first occurrence of `text[pos-min_len..pos]`,
    /// but this is *not* the first occurrence of `text[pos-min_len+1..pos]`,
    /// i.e., `s[1..]A` has a previous occurrence.
    min_len: usize,
    /// The length of the longest `sA` such that `sA` is an RME.
    /// Adding one extra character on the left will already be 'anchored' elsewhere.
    /// That is, `text[pos-max_len..pos-1]` is not leftmost here, but
    /// `text[pos-max_len-1..pos-1]` (and longer) _are_ leftmost here.
    max_len: usize,
    /// The position in the text (before `pos`) where the leftmost occurrence of `text[pos-min_len+1..pos]` ends.
    /// Note this this position may or may not corresponds to an RME.
    suffix_pos: usize,
    /// The position in the text (before `suffix_pos`) where the leftmost occurrence of `text[pos-min_len+1..pos]` is anchored.
    suffix_anchor_pos: usize,
}

impl Stpd {
    pub fn new(text: &[u8]) -> Self {
        // 1. Build SA
        // 2. Iterate text left to right, and mark samples as needed.
        // 3.
        todo!();
    }

    /// Find the range of `spa` that has `q` as a suffix.
    // TODO: Do everything in reverse instead, so we have normal lex comparisons?
    fn binary_search(&self, q: &[u8]) -> Range<usize> {
        let start = self
            .spa
            .binary_search_by(|rme| match cmp_colex(&self.text[..rme.pos], q) {
                Ordering::Equal => Ordering::Greater,
                x => x,
            })
            .unwrap_err();
        let end = self
            .spa
            .binary_search_by(|rme| match cmp_colex(&self.text[..rme.pos], q) {
                Ordering::Equal => Ordering::Less,
                x => x,
            })
            .unwrap_err();
        start..end
    }

    /// Return the leftmost sampled RME matching q, if it exists.
    fn search_anchor(&self, q: &[u8]) -> Option<&Anchor> {
        assert!(!q.is_empty());
        let range = self.binary_search(q);
        if range.is_empty() {
            return None;
        }
        // Find the smallest index in the range.
        // TODO: O(1) RMQ?
        let rme_index = range.start
            + self.spa[range]
                .iter()
                .position_min_by_key(|rme| rme.pos)
                .unwrap();
        Some(&self.spa[rme_index])
    }

    /// Leftmost occurrence of q, if it occurs.
    pub fn locate_one(&self, q: &[u8]) -> Option<Range<usize>> {
        let m = self.extend(q, 0..0, &self.spa[0]).0;
        if m.len() == q.len() {
            Some(m)
        } else {
            None
        }
    }

    /// Given q, and the location of an already matched prefix, find the longest prefix that occurs and return:
    /// - the leftmost slice of the text corresponding to the matched prefix,
    /// - the last RME anchor.
    pub fn extend<'s>(
        &'s self,
        q: &[u8],
        prefix_match: Range<usize>,
        mut anchor: &'s Anchor,
    ) -> (Range<usize>, &'s Anchor) {
        assert_eq!(self.text[prefix_match.clone()], q[..prefix_match.len()]);
        assert!(prefix_match.start <= anchor.pos && anchor.pos <= prefix_match.end);

        // Number of matched characters of q.
        let mut i = prefix_match.len();
        // *End* position of leftmost occurrence in text of q[..i].
        let mut pos = prefix_match.end;
        while i < q.len() {
            if self.text[pos] == q[i] {
                i += 1;
                pos += 1;
                continue;
            }

            // q[..i] does not occur at `pos`, so is either an RME or does not occur at all.
            let Some(new_anchor) = self.search_anchor(&q[..=i]) else {
                break;
            };
            anchor = new_anchor;
            pos = anchor.pos;
            i += 1;
        }
        (pos - i..pos, anchor)
    }

    /// Given the leftmost occurrence ending at `pos` (inclusive!) with the given RME anchor
    /// (which might end before pos), find the suffixlink, ie, the longest
    /// suffix of `text[..pos]` that occurs earlier in the text.
    ///
    ///                        anchor/RME
    /// extra       minlen     pos        pos   unmatched
    /// ...Z   AB   CDEFGH     I          JK    X...
    ///        ----------------------------- current matched text
    ///              ----------------------- suffix
    ///              -----------
    ///
    ///     suffix anchor RME
    ///     v
    /// DEF G HI JY
    ///        ^ suffix pos
    ///
    /// We want to find occurrence that also includes X.
    /// But extending here with a binary search for ...X does not work,
    /// so ABCDEFGHIJKX does not occur.
    /// But we do already get at lower bound on the length of the final match.
    /// `minlen` tells us that this is the leftmost CDEFGHI, so the
    /// suffix link of dropping A stays at the same anchor, as does dropping B
    /// after that. Dropping the C *does* give a new leftmost occurrence of the anchor string DEFGHI.
    ///
    /// `suffix_pos` is the position of the final I of the leftmost occurrence of DEFGHI.
    /// The corresponding anchor might be somewhat to the left of it though.
    /// From there, do a normal while-loop of extend and binary search to try to match JKX.
    /// This might not work, in which case take the suffix link of the last anchor and repeat.
    ///         

    /// `m`: The currently matched substring, of which we want the first non-trivial suffix link.
    /// `anchor`: The corresponding anchor.
    fn suffix_link<'s>(
        &'s self,
        mut matched: Range<usize>,
        mut anchor: &'s Anchor,
    ) -> (&'s Anchor, usize) {
        let mut target = &self.text[matched.clone()];

        loop {
            // anchor itself is in the range
            assert!(matched.start <= anchor.pos && anchor.pos <= matched.end);
            // m has the right length to be anchored here
            assert!(
                anchor.pos - anchor.max_len <= matched.start
                    && matched.start <= anchor.pos - anchor.min_len
            );
            // Shrink the target to the next suffix-link length.
            target = &target[(anchor.pos - anchor.min_len + 1) - matched.start..];
            // The range of text matched by the suffix link.
            matched = anchor.suffix_pos - anchor.min_len + 1..anchor.suffix_pos;
            // The suffix link anchor.
            anchor = self
                .search_anchor(&self.text[matched.start..anchor.suffix_anchor_pos])
                .unwrap();
            // Extend the match as much as possible.
            (matched, anchor) = self.extend(target, matched, anchor);
            assert_eq!(&self.text[matched.clone()], &target[..matched.len()]);
            // If the entire target suffix matched, return it.
            if matched.len() == target.len() {
                return (anchor, anchor.pos);
            }
            // Otherwise, shrink the target further.
        }
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
