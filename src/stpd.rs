use std::{
    bstr::{ByteStr, ByteString},
    cmp::{Ordering, Reverse},
    ops::Range,
};

/// Implementation of lex-min STPD.
// TODO: Suffix lookup for binary search.
pub struct Stpd {
    /// The input text.
    // TODO: Relative LZ encoding.
    text: ByteString,
    /// STPD samples, covering all left-most occurrences of all right-maximal extensions.
    /// Sorted in co-lex order.
    /// spa: sampled prefix array
    /// `spa[0]` is the 'trivial' empty anchor at the root of the text.
    // TODO: BTree instead so we can insert samples.
    // TODO: Efficient range min.
    spa: tiered_vector::Vector<Anchor>,

    /// The longest suffix `text[pos-seen_before.len()..pos]` that was seen before at `text[seen_before]`.
    seen_before: Range<usize>,
    /// The anchor index of `text[seen_before]`.
    anchor_idx: usize,
}

/// A right maximal extension is a substring sA of T such that s occurs at least twice in the text.
/// We do path compression, which means that the A only sampled if `sA` is
/// leftmost, but `s` has another occurrence further to the left.
/// `pos` indicates the length of the leftmost text prefix `text[..pos]` ending in `sA`.
/// This position is an RME for an interval of length of `s`:
/// - For very short `s'`, `s'A` will already occur to the left of `pos`.
/// - For very long `s''`, this will already be the leftmost occurrence of `s''`.
// TODO: Inline the last u64 of the prefix.
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
    // /// The length of the longest `sA` such that `sA` is an RME.
    // /// Adding one extra character on the left will already be 'anchored' elsewhere.
    // /// That is, `text[pos-max_len..pos-1]` is not leftmost here, but
    // /// `text[pos-max_len-1..pos-1]` (and longer) _are_ leftmost here.
    // max_len: usize,
    /// The position in the text (before `pos`) where the leftmost occurrence of `text[pos-min_len+1..pos]` ends.
    /// Note this this position may or may not corresponds to an RME.
    suffix_pos: usize,
    /// The position in the text (before `suffix_pos`) where the leftmost occurrence of `text[pos-min_len+1..pos]` is anchored.
    /// *Not* the index in `STPD::spa` of the anchor itself, since this can change over time.
    suffix_anchor_pos: usize,
}

impl Stpd {
    pub fn new(text: &[u8]) -> Self {
        log::debug!("NEW");
        let mut stpd = Self {
            text: ByteString(vec![]),
            spa: tiered_vector::Vector::new(),
            seen_before: 0..0,
            anchor_idx: 0,
        };
        stpd.spa.push(Anchor {
            pos: 0,
            min_len: 0,
            // max_len: 0,
            suffix_pos: 0,
            suffix_anchor_pos: 0,
        });
        stpd.push(text);
        stpd
    }
    pub fn push(&mut self, text: &[u8]) {
        log::info!("Push {}", ByteStr::new(text));
        let old_len = self.text.len();
        self.text.extend_from_slice(text);
        let text = &self.text;

        let mut seen_before = self.seen_before.clone();
        let mut anchor_idx = self.anchor_idx;
        // No need to create an anchor for the first character.
        for pos in old_len.max(1)..text.len() {
            // Append text[pos].
            let c = text[pos];
            debug_assert_eq!(
                text[seen_before.clone()],
                text[pos - seen_before.len()..pos]
            );

            let extended = &text[pos - seen_before.len()..=pos];
            log::debug!("Extended: {extended:?}");

            // the seen-before part with one extra character

            // Found by extending match at current anchor.
            if text[seen_before.end] == c {
                log::debug!("Found by extending current match.");
                seen_before.end += 1;
                continue;
            }
            log::info!(
                "Pos {pos} Push {}. Seen before: {}=|{seen_before:?}| with anchor {}",
                c as char,
                seen_before.len(),
                self.spa[anchor_idx].pos
            );

            // Search prefix array for match with additional character.
            if let Some(((ai, _anchor), sb)) = self.search_anchor(extended) {
                log::debug!("Found match of {extended:?} via binary search at {sb:?}.");
                anchor_idx = ai;
                seen_before = sb;
                log::debug!("New anchor {anchor_idx}");
                debug_assert_eq!(
                    &text[seen_before.clone()],
                    &text[pos - seen_before.len() + 1..=pos]
                );
                continue;
            }

            log::debug!("Add new anchor for max_len={}", extended.len());

            // This is the first occurrence of `extended`, so we add an anchor for it.
            // We now search for the "min_len", ie the longest suffix of
            // `extended` that occurs before, and its anchor.
            // let max_len = extended.len();
            let anchor;
            ((anchor_idx, anchor), seen_before) =
                self.longest_existing_suffix(seen_before, anchor_idx, pos);

            let new_anchor = Anchor {
                pos: pos + 1,
                // max_len,
                min_len: seen_before.len() + 1,
                suffix_pos: seen_before.end,
                suffix_anchor_pos: anchor.pos,
            };

            log::info!("New anchor {new_anchor:?}");

            // Insert anchor
            let pos_range = self.binary_search(&text[..=pos]);
            assert_eq!(pos_range.start, pos_range.end);
            let insert_idx = pos_range.start;
            log::debug!(
                "Insert anchor #{} at position {}",
                self.spa.len(),
                insert_idx
            );
            self.spa.insert(insert_idx, new_anchor);
            if insert_idx <= anchor_idx {
                log::debug!(
                    "Inserted before current anchor idx; updating that to {}",
                    anchor_idx + 1
                );
                anchor_idx += 1;
            }
        }
        log::error!("Num anchors: {}", self.spa.len());
        self.seen_before = seen_before;
        self.anchor_idx = anchor_idx;
    }

    fn longest_existing_suffix(
        &self,
        seen_before: Range<usize>,
        anchor_idx: usize,
        pos: usize,
    ) -> ((usize, &Anchor), Range<usize>) {
        #[cfg(not(debug_assertions))]
        if seen_before.len() > 20 {
            // self.longest_existing_suffix_via_binary_search(seen_before, anchor_idx, pos)
            self.longest_existing_suffix_via_exponential_search(seen_before, anchor_idx, pos)
        } else {
            self.longest_existing_suffix_via_links(seen_before, anchor_idx, pos)
        }

        #[cfg(debug_assertions)]
        {
            let a = self.longest_existing_suffix_via_binary_search(
                seen_before.clone(),
                anchor_idx,
                pos,
            );
            let b = self.longest_existing_suffix_via_links(seen_before.clone(), anchor_idx, pos);
            let c = self.longest_existing_suffix_via_exponential_search(
                seen_before.clone(),
                anchor_idx,
                pos,
            );
            assert_eq!(a.0 .0, b.0 .0);
            assert_eq!(a.0 .0, c.0 .0);
            assert_eq!(a.1, b.1);
            assert_eq!(a.1, c.1);
            a
        }
    }

    /// Given that the current suffix occurs before at `seen_before` with the given `anchor_idx`,
    /// find the longest suffix of `text[..=pos]` that occurs before, and its anchor.
    ///
    /// This version binary searches the length of the suffix to find the longest that occurs.
    fn longest_existing_suffix_via_binary_search(
        &self,
        seen_before: Range<usize>,
        _anchor_idx: usize,
        pos: usize,
    ) -> ((usize, &Anchor), Range<usize>) {
        let extended = &self.text[pos - seen_before.len()..=pos];

        let mut l = 0;
        let mut h = extended.len();
        log::warn!("Binary search for suffix len 0..{h}");
        while l < h {
            let mid = (l + h + 1) / 2;
            log::debug!("l {l} mid {mid} h {h}");
            let q = &extended[extended.len() - mid..];

            if let Some(mut m) = self.locate_one(q) {
                assert_eq!(m.len(), mid);
                // Try to extend the match on the left to increase l.
                debug_assert_eq!(&self.text[m.clone()], q);
                l = mid;
                while l + 1 < extended.len() && m.start > 0 {
                    let c = extended[extended.len() - 1 - l];
                    if self.text[m.start - 1] != c {
                        break;
                    }
                    l += 1;
                    m.start -= 1;
                }
                assert_eq!(m.len(), l);
                // assert_eq!(&self.text[m], &extended[extended.len() - l..]);
                if l > mid {
                    log::info!("Grow l from {mid} to {l}");
                }
            } else {
                h = mid - 1;
            }
        }
        log::info!("l {l} h {h}");
        let q = &extended[extended.len() - l..];
        let (seen_before, (anchor_idx, anchor)) = self.extend(q, 0..0, 0, &self.spa[0]);
        assert_eq!(seen_before.len(), q.len());
        ((anchor_idx, anchor), seen_before)
    }

    /// Exponential search version that doubles the prefix length.
    fn longest_existing_suffix_via_exponential_search(
        &self,
        seen_before: Range<usize>,
        _anchor_idx: usize,
        pos: usize,
    ) -> ((usize, &Anchor), Range<usize>) {
        let extended = &self.text[pos - seen_before.len()..=pos];

        let mut l = 0;
        let mut h = extended.len();
        log::warn!("exponential search for suffix len 0..{h}");
        while l < h {
            let mid = if 2 * l + 1 < h {
                if h < 20 {
                    2 * l + 1
                } else {
                    (2 * l + 1).max(7)
                }
            } else {
                (l + h + 1) / 2
            };
            log::debug!("l {l} mid {mid} h {h}");
            let q = &extended[extended.len() - mid..];

            if let Some(mut m) = self.locate_one(q) {
                assert_eq!(m.len(), mid);
                // Try to extend the match on the left to increase l.
                debug_assert_eq!(&self.text[m.clone()], q);
                l = mid;
                while l + 1 < extended.len() && m.start > 0 {
                    let c = extended[extended.len() - 1 - l];
                    if self.text[m.start - 1] != c {
                        break;
                    }
                    l += 1;
                    m.start -= 1;
                }
                assert_eq!(m.len(), l);
                // assert_eq!(&self.text[m], &extended[extended.len() - l..]);
                if l > mid {
                    log::info!("Grow l from {mid} to {l}");
                }
            } else {
                h = mid - 1;
            }
        }
        log::debug!("l {l} h {h}");
        let q = &extended[extended.len() - l..];
        let (seen_before, (anchor_idx, anchor)) = self.extend(q, 0..0, 0, &self.spa[0]);
        assert_eq!(seen_before.len(), q.len());
        ((anchor_idx, anchor), seen_before)
    }

    /// Given that the current suffix occurs before at `seen_before` with the given `anchor_idx`,
    /// find the longest suffix of `text[..=pos]` that occurs before, and its anchor.
    ///
    /// This version repeatedly follows suffix links.
    fn longest_existing_suffix_via_links(
        &self,
        mut seen_before: Range<usize>,
        mut anchor_idx: usize,
        pos: usize,
    ) -> ((usize, &Anchor), Range<usize>) {
        let c = self.text[pos];
        let extended = &self.text[pos - seen_before.len()..=pos];

        // Suppose the existing longest seen before suffix is ABCDEF and it matches
        // around some anchor as ABC.DEF.
        // We want to find whether there is a previous occurrence of ABCDEFG, and if not, the longest suffix for which there is.
        // 3 options for the leftmost longest suffix match:
        // 1) at the current anchor. => The new suffix is not an RME and already handled above.
        // 2) right of the anchor. The G will be sampled as an RME, and we find it via binary search.
        // 3) left of the anchor. We find it by repeated suffix-link and extend.

        // 2) Find the longest match in the prefix array.
        let range = self.binary_search(&self.text[..=pos]);
        assert!(range.is_empty());
        // Test the strings before and after.
        // Find all that have the maximal LCS length, and of those, report the leftmost.
        let mut max_lcs = (0, Reverse(usize::MAX), 0);
        // (lcs length, text index, anchor idx)

        let mut i = range.start;
        while i > 0 {
            i -= 1;
            let lcs_len = lcs(&self.text[..self.spa[i].pos], extended);
            if lcs_len < max_lcs.0 {
                break;
            }
            max_lcs = max_lcs.max((lcs_len, Reverse(self.spa[i].pos), i));
        }
        i = range.start;
        while i < self.spa.len() {
            let lcs_len = lcs(&self.text[..self.spa[i].pos], extended);
            if lcs_len < max_lcs.0 {
                break;
            }
            max_lcs = max_lcs.max((lcs_len, Reverse(self.spa[i].pos), i));
            i += 1;
        }
        let right_anchor_idx = max_lcs.2;
        let right_anchor = &self.spa[right_anchor_idx];
        let right_seen_before = right_anchor.pos - max_lcs.0..right_anchor.pos;
        log::info!(
            "Right seen before: {right_seen_before:?} with anchor {}",
            right_anchor.pos
        );

        // 3) Repeatedly take suffix links to find the min_len, ie the length of the shortest suffix
        // for which the just pushed character is the anchor.
        let mut anchor = &self.spa[anchor_idx];
        // TODO: Prune search once it's worse than what we see on the right.
        log::warn!(
            "({pos}) Suffix link of {}=|{seen_before:?}| extended by {}",
            seen_before.len(),
            c as char
        );
        while seen_before.len() > 0 {
            // Seen before is one less than that.
            // Take suffix link of the anchor.
            ((anchor_idx, anchor), seen_before) = self.suffix_link(seen_before.clone(), anchor);
            debug_assert_eq!(
                self.text[seen_before.clone()],
                self.text[pos - seen_before.len()..pos]
            );

            // If we can extend the suffix match with the right character, we update `seen_before` and are done.
            if self.text[seen_before.end] == c {
                log::debug!("Found previous occurrence by extending suffix link.");
                seen_before.end += 1;
                break;
            }

            // Otherwise, we might be able to find an existing RME anchor.
            let suffix = &self.text[pos - seen_before.len()..=pos];
            log::debug!("Binary search for suffix {suffix:?}",);
            if let Some(((ai, a), sb)) = self.search_anchor(suffix) {
                log::debug!("Found previous occurrence via binary search at {sb:?} {a:?}.");
                anchor_idx = ai;
                anchor = a;
                seen_before = sb;
                debug_assert_eq!(
                    &self.text[seen_before.clone()],
                    &self.text[pos + 1 - seen_before.len()..=pos]
                );
                assert!(seen_before.len() <= right_seen_before.len());
                break;
            }
            log::debug!("Take another suffix link");

            // Otherwise, this suffix cannot be extended with `c`, and we take further suffix links.
        }

        // Take the max of the two options.
        if right_seen_before.len() > seen_before.len()
            || (right_seen_before.len() == seen_before.len() && right_anchor.pos < anchor.pos)
        {
            log::debug!("Found better match on the right.");
            anchor = right_anchor;
            seen_before = right_seen_before;
            anchor_idx = right_anchor_idx;
        }
        ((anchor_idx, anchor), seen_before)
    }

    /// `cmp`: true when anchor < query.
    fn binary_search_by(&self, cmp: impl Fn(&Anchor) -> bool) -> usize {
        let mut l = 0;
        let mut h = self.spa.len();
        while l < h {
            let m = (l + h) / 2;
            if cmp(&self.spa[m]) {
                l = m + 1;
            } else {
                h = m;
            }
        }
        l
    }

    /// Find the range of `spa` that has `q` as a suffix.
    fn binary_search(&self, q: &[u8]) -> Range<usize> {
        let start =
            self.binary_search_by(|rme| cmp_colex(&self.text[..rme.pos], q) == Ordering::Less);
        let end =
            self.binary_search_by(|rme| cmp_colex(&self.text[..rme.pos], q) != Ordering::Greater);
        start..end
    }

    /// Return the leftmost sampled RME matching q and the matching range, if it exists.
    fn search_anchor(&self, q: &[u8]) -> Option<((usize, &Anchor), Range<usize>)> {
        let range = self.binary_search(q);
        if range.is_empty() {
            return None;
        }
        // Find the smallest index in the range.
        // TODO: O(1) RMQ?
        let rme_index = range
            .into_iter()
            .min_by_key(|idx| self.spa[*idx].pos)
            .unwrap();
        let anchor = &self.spa[rme_index];
        let range = anchor.pos - q.len()..anchor.pos;
        debug_assert_eq!(&self.text[range.clone()], q);
        Some(((rme_index, anchor), range))
    }

    /// Leftmost occurrence of q, if it occurs.
    pub fn locate_one(&self, q: &[u8]) -> Option<Range<usize>> {
        let m = self.extend(q, 0..0, 0, &self.spa[0]).0;
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
        mut anchor_idx: usize,
        mut anchor: &'s Anchor,
    ) -> (Range<usize>, (usize, &'s Anchor)) {
        debug_assert_eq!(self.text[prefix_match.clone()], q[..prefix_match.len()]);
        assert!(prefix_match.start <= anchor.pos && anchor.pos <= prefix_match.end);

        // Number of matched characters of q.
        let mut i = prefix_match.len();
        // *End* position of leftmost occurrence in text of q[..i].
        let mut pos = prefix_match.end;
        let mut searches = 0;
        while i < q.len() {
            if unsafe { *self.text.get_unchecked(pos) } == q[i] {
                i += 1;
                pos += 1;
                continue;
            }

            searches += 1;

            // q[..i] does not occur at `pos`, so is either an RME or does not occur at all.
            let Some(((new_anchor_idx, new_anchor), _)) = self.search_anchor(&q[..=i]) else {
                break;
            };
            anchor_idx = new_anchor_idx;
            anchor = new_anchor;
            pos = anchor.pos;
            i += 1;
        }

        let range = pos - i..pos;
        log::info!(
            "extend |q|={} with {searches} searches from {}=|{prefix_match:?}| to {}=|{range:?}| anchored at {}",
            q.len(),
            prefix_match.len(),
            range.len(),
            anchor.pos
        );
        (range, (anchor_idx, anchor))
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
    ) -> ((usize, &'s Anchor), Range<usize>) {
        assert!(matched.len() > 0);
        let mut target = &self.text[matched.clone()];
        log::info!("Suffix link of: matched={matched:?} {anchor:?} ");

        let mut anchor_idx;

        loop {
            // anchor itself is in the range
            assert!(matched.start <= anchor.pos && anchor.pos <= matched.end);
            // m has the right length to be anchored here
            assert!(
                // anchor.pos - anchor.max_len <= matched.start &&
                matched.start <= anchor.pos - anchor.min_len
            );
            // Special case: suffix links of the root anchor go to itself.
            if anchor.pos == 0 {
                // Shrink the target to the next suffix-link length.
                target = &target[1..];
                // The range of text matched by the suffix link.
                matched = 0..0;
            } else {
                // Shrink the target to the next suffix-link length.
                target = &target[(anchor.pos + 1 - (anchor.min_len - 0)) - matched.start..];
                // The range of text matched by the suffix link.
                matched = anchor.suffix_pos + 1 - (anchor.min_len - 0)..anchor.suffix_pos;
            }
            log::debug!("Shrink target to {target:?}");
            // log::info!("Updated matched to {matched:?}");
            // The suffix link anchor.
            (anchor_idx, anchor) = self
                .search_anchor(&self.text[matched.start..anchor.suffix_anchor_pos])
                .unwrap()
                .0;
            log::debug!("suffix link anchor {anchor:?}");
            // Extend the match as much as possible.
            log::debug!("Extend {target:?} {matched:?} {anchor:?}");
            (matched, (anchor_idx, anchor)) = self.extend(target, matched, anchor_idx, anchor);
            log::debug!("Extend into {matched:?} {anchor:?}");
            debug_assert_eq!(&self.text[matched.clone()], &target[..matched.len()]);
            // If the entire target suffix matched, return it.
            if matched.len() == target.len() {
                log::debug!("Found suffix link: matched={matched:?} anchor={anchor:?}");
                return ((anchor_idx, anchor), matched);
            }
            log::debug!(
                "Current match {matched:?} = {} is not yet full target {target:?}",
                &self.text[matched.clone()]
            );
            // Otherwise, shrink the target further.
        }
        // return (anchor, matched);
    }
}

/// Length of longest common suffix.
fn lcs(a: &[u8], b: &[u8]) -> usize {
    // TODO: u64-based comparisons.
    let mut i = 0;
    let min = Ord::min(a.len(), b.len());
    while i < min && a[a.len() - 1 - i] == b[b.len() - 1 - i] {
        i += 1;
    }
    return i;
}

/// Read u64 ending at position i >= 8.
fn read_u64(text: &[u8], i: usize) -> u64 {
    debug_assert!(i >= 8);
    unsafe { u64::from_le_bytes(text.get_unchecked(i - 8..i).try_into().unwrap()) }
}

fn read_last_u64(text: &[u8], i: usize, len: usize) -> u64 {
    debug_assert!(len < 8);
    let mut data = [0; 8];
    data[8 - len..].copy_from_slice(unsafe { text.get_unchecked(i - len..i) });
    u64::from_le_bytes(data)
}

/// co-lex compare q with a text prefix.
/// Returns `Equal` when `q` is a suffix of `text`.
fn cmp_colex(text: &[u8], q: &[u8]) -> Ordering {
    let min_len = Ord::min(text.len(), q.len());
    let min = min_len / 8 * 8;
    for i in (0..min).step_by(8) {
        let text_chunk = read_u64(text, text.len() - i);
        let q_chunk = read_u64(q, q.len() - i);
        if text_chunk != q_chunk {
            return Ord::cmp(&text_chunk, &q_chunk);
        }
    }
    let text_chunk = read_last_u64(text, text.len() - min, min_len - min);
    let q_chunk = read_last_u64(q, q.len() - min, min_len - min);
    if text_chunk != q_chunk {
        return Ord::cmp(&text_chunk, &q_chunk);
    }
    if min_len == q.len() {
        return Ordering::Equal;
    }
    if min_len == text.len() {
        return Ordering::Less;
    }
    unreachable!();
}
