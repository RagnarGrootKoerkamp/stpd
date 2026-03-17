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
#[derive(Debug)]
pub struct Anchor {
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
    pub fn new(text: Vec<u8>) -> Self {
        // 1. Build SA
        // 2. Iterate text left to right, and mark samples as needed.
        // 3.

        let mut stpd = Self {
            text,
            spa: vec![Anchor {
                pos: 0,
                min_len: 0,
                max_len: 0,
                suffix_pos: usize::MAX,
                suffix_anchor_pos: usize::MAX,
            }],
        };
        log::info!("New");

        let text = &stpd.text;

        // The substring `text[pos-seen_before.len()..pos]` was seen before at `text[seen_before]`.
        // Also, `text[pos-seen_before.len()-1..pos]` was *not* seen before.
        // `seen_before.end < pos`.
        let mut seen_before = 0..0;
        // The anchor of `text[seen_before]`.
        let mut anchor_idx = 0;

        for pos in 1..text.len() {
            // Append text[pos].
            let c = text[pos];
            log::info!(
                "Push {c}. Seen before: text[{seen_before:?}]={:?} with anchor {anchor_idx}",
                &text[seen_before.clone()]
            );

            let extended = &text[pos - seen_before.len()..=pos];
            log::info!("Extended: {extended:?}");

            // the seen-before part with one extra character

            // Found by extending match at current anchor.
            if text[seen_before.end] == c {
                log::info!("Found by extending current match.");
                seen_before.end += 1;
                continue;
            }

            // Search prefix array for match with additional character.
            if let Some((anchor, sb)) = stpd.search_anchor(extended) {
                log::info!("Found via binary search.");
                anchor_idx = stpd.spa.element_offset(anchor).unwrap();
                seen_before = sb;
                log::info!("New anchor {anchor_idx} with seen_before={seen_before:?}");
                assert_eq!(
                    &text[seen_before.clone()],
                    &text[pos - seen_before.len()..pos]
                );
                continue;
            }

            log::info!("Add new anchor for max_len={}", extended.len());

            // This is the first occurrence of `extended`, so we add an anchor for it.
            let mut new_anchor = Anchor {
                pos,
                max_len: extended.len(),
                // TODO
                min_len: usize::MAX,
                suffix_pos: usize::MAX,
                suffix_anchor_pos: usize::MAX,
            };

            // TODO: Repeatedly take suffix links to find the min_len, ie the length of the shortest suffix
            // for which the just pushed character is the anchor.
            let mut anchor = &stpd.spa[anchor_idx];
            loop {
                // Seen before is one less than that.
                // Take suffix link of the anchor.
                (anchor, seen_before) = stpd.suffix_link(seen_before, anchor);
                anchor_idx = stpd.spa.element_offset(anchor).unwrap();

                // If we can extend the suffix match with the right character, we update `seen_before` and are done.
                if text[seen_before.end] == c {
                    log::info!("Found previous occurrence by extending suffix link.");
                    seen_before.end += 1;
                    break;
                }

                // Otherwise, we might be able to find an existing RME anchor.
                let suffix = &text[pos - seen_before.len()..=pos];
                if let Some((a, sb)) = stpd.search_anchor(suffix) {
                    log::info!("Found previous occurrence via binary search.");
                    anchor = a;
                    seen_before = sb;
                    anchor_idx = stpd.spa.element_offset(anchor).unwrap();
                    assert_eq!(
                        &text[seen_before.clone()],
                        &text[pos - seen_before.len()..pos]
                    );
                    break;
                }
                log::info!("Take another suffix link");

                // Otherwise, this suffix cannot be extended with `c`, and we take further suffix links.
            }

            // Update anchor
            new_anchor.min_len = seen_before.len() + 1;
            new_anchor.suffix_pos = seen_before.end;
            new_anchor.suffix_anchor_pos = anchor.pos;

            log::info!("New anchor {new_anchor:?}");

            // Insert anchor
            let pos_range = stpd.binary_search(&text[..pos]);
            assert_eq!(pos_range.start, pos_range.end);
            log::info!("Insert new anchor at position {}", pos_range.start);
            stpd.spa.insert(pos_range.start, new_anchor);
        }
        stpd
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

    /// Return the leftmost sampled RME matching q and the matching range, if it exists.
    fn search_anchor(&self, q: &[u8]) -> Option<(&Anchor, Range<usize>)> {
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
        let anchor = &self.spa[rme_index];
        let range = anchor.pos - q.len()..anchor.pos;
        Some((anchor, range))
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
            let Some((new_anchor, _)) = self.search_anchor(&q[..=i]) else {
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
    /// Returns:
    /// - the index in `self.spa` of the anchor
    /// - the matched range of text, which equals a suffix of `matched`.
    fn suffix_link<'s>(
        &'s self,
        mut matched: Range<usize>,
        mut anchor: &'s Anchor,
    ) -> (&'s Anchor, Range<usize>) {
        let mut target = &self.text[matched.clone()];
        log::info!("Suffix link of: matched={matched:?} anchor={anchor:?} target={target:?}");

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
                .unwrap()
                .0;
            // Extend the match as much as possible.
            (matched, anchor) = self.extend(target, matched, anchor);
            assert_eq!(&self.text[matched.clone()], &target[..matched.len()]);
            // If the entire target suffix matched, return it.
            if matched.len() == target.len() {
                log::info!("Found suffix link: matched={matched:?} anchor={anchor:?}");
                return (anchor, matched);
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
