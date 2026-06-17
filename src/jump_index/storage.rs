#![allow(unused_imports)]
use std::{
    cmp::max,
    collections::{BTreeSet, HashMap},
};

use itertools::{Either, Itertools};
use mem_dbg::MemSize;
#[cfg(feature = "mphf")]
use ptr_hash::{hash::FastIntHash, DefaultPtrHash, PtrHash, PtrHashParams};
use sux::{
    bits::BitVec,
    dict::EliasFanoBuilder,
    rank_sel::SelectAdaptConst,
    traits::{Pred, Succ},
};
use voracious_radix_sort::RadixSort;

use crate::gbs;

use super::{
    build_ef,
    link::{BareEf, Link, LinkEf, SuffixLink},
};

#[cfg(feature = "mphf")]
#[derive(MemSize)]
pub struct MphfStore {
    /// PHF mapping (source, c) to index.
    phf: DefaultPtrHash<FastIntHash, usize>,
    /// Bitvec indicating for each slot.
    /// 1 per slot
    /// 0 before the 1 for every additional link target here.
    bits: SelectAdaptConst<BitVec>,
    /// targets
    targets: Vec<u32>,
}

#[cfg(feature = "mphf")]
impl MphfStore {
    pub fn new(links: &LinkEf) -> Self {
        let keys: Vec<_> = links
            .iter()
            .map(|x| Link::from_key(x).source_c())
            .dedup()
            .collect();
        let mut params = PtrHashParams::default();
        params.remap = false;
        params.alpha = 0.98;
        eprintln!("Building PHF..");
        let phf = DefaultPtrHash::new(&keys, params);
        let slots = phf.max_index();
        let mut lens = vec![0u8; slots];
        for link in links.iter() {
            let link = Link::from_key(link);
            let idx = phf.index_no_remap(&link.source_c());
            lens[idx] += 1;
        }
        // increment all lens to at least 1
        for len in &mut lens {
            *len = 1.max(*len);
        }
        let target_slots = lens.iter().map(|x| *x as usize).sum();

        // build bitvec. 1 + 0*(len-1) for each slot, + 1 at the end.
        eprintln!("Build bitvec..");
        let mut bits = BitVec::with_capacity(lens.iter().map(|x| *x as usize).sum());
        for len in &lens {
            bits.push(true);
            bits.extend(std::iter::repeat(false).take((*len - 1) as usize));
        }
        bits.push(true);
        assert_eq!(bits.len(), target_slots + 1);
        let bits = SelectAdaptConst::<_, _>::new(bits);

        eprintln!("Build targets..");
        lens.fill(0);
        let mut targets = vec![0u32; target_slots];
        for link in links.iter() {
            let link = Link::from_key(link);
            let idx = phf.index_no_remap(&link.source_c());
            let pos = lens[idx] as usize;
            targets[pos] = link.target() as u32;
            lens[idx] += 1;
        }
        MphfStore { phf, bits, targets }
    }

    pub fn size(&self) -> String {
        format!(
            "MphfStore:\n  phf: {:.3} GB\n  bits: {:.3} GB\n  targets: {:.3} GB",
            gbs(&self.phf),
            gbs(&self.bits),
            gbs(&self.targets)
        )
    }
}

#[derive(MemSize)]
pub struct RelativeStore {
    fwd_ef: LinkEf,
    suf_ef: LinkEf,
    /// (source, c) -> (len, index) with fwd_targets[len][len*index .. len*(index+1)]
    fwd_idx: HashMap<usize, (u32, u32)>,
    /// Index in corresponding fwd_targets that are free.
    free_targets: Vec<Vec<u32>>,
    /// lcp and target
    fwd_targets: Vec<Vec<(u32, u32)>>,
    suf_set: Vec<u128>,
}

impl RelativeStore {
    pub fn new(fwd_ef: LinkEf, suf_ef: LinkEf) -> Self {
        Self {
            fwd_ef,
            suf_ef,
            fwd_idx: HashMap::new(),
            free_targets: vec![vec![]; 40],
            fwd_targets: vec![vec![]; 40],
            suf_set: vec![],
        }
    }
    pub fn fwd_len(&self) -> usize {
        // FIXME: Account for free slots.
        self.fwd_ef.len() + self.fwd_targets.iter().map(|x| x.len()).sum::<usize>()
    }
    pub fn suf_len(&self) -> usize {
        self.suf_ef.len() + self.suf_set.len()
    }
    pub fn insert_fwd(&mut self, link: Link) {
        // eprintln!("INSERT {link:?}");
        // get current len
        // if new, append to len 1
        // otherwise, move data to len+1
        // and fill gap left behind in previous data
        let source_c = link.source_c();
        let val = (link.lcp() as u32, link.target() as u32);
        let entry = self.fwd_idx.entry(source_c).or_insert((0, 0));
        let old_len = entry.0 as usize;
        entry.0 += 1;
        let len = entry.0 as usize;

        // make space in next layer
        let idx = if let Some(idx) = self.free_targets[len as usize].pop() {
            idx as usize
        } else {
            self.fwd_targets[len as usize].extend(std::iter::repeat((0, 0)).take(len as usize));
            self.fwd_targets[len as usize].len() - len
        };
        let old_idx = entry.1 as usize;
        entry.1 = idx as u32;

        // copy from existing layer
        for i in 0..(len - 1) {
            let old_val = self.fwd_targets[old_len][old_idx + i];
            self.fwd_targets[len][idx + i] = old_val;
        }
        // push on next layer
        self.fwd_targets[len][idx + len - 1] = val;
        // sort values
        self.fwd_targets[len][idx..idx + len].sort_unstable_by_key(|x| x.0);

        // Mark previous slot as free.
        if old_len > 0 {
            self.free_targets[old_len].push(old_idx as u32);
        }
    }
    pub fn insert_suf(&mut self, link: SuffixLink) {
        let key = link.key();
        if self.suf_set.is_empty() {
            let last_in_ef = self.suf_ef.iter_back().next().unwrap_or(0);
            assert!(
                key >= last_in_ef,
                "Link {link:?}={key} is not larger than last {:?}={last_in_ef}",
                SuffixLink::from_key(last_in_ef)
            );
        } else {
            assert!(key > *self.suf_set.last().unwrap());
        }
        self.suf_set.push(key);
    }
    /// Return the first fwd link >= `link`, but only if it has matching source and character.
    pub fn fwd_succ(&self, link: Link) -> Option<Link> {
        // eprintln!("SUCC OF {link:?}");
        let x = link.key();
        let k1 = self.fwd_ef.succ(x).map(|x| Link::from_key(x.1));
        let k2 = try {
            let &(len, idx) = self.fwd_idx .get(&link.source_c())?;
            let len = len as usize;
            let idx = idx as usize;
            assert!(len > 0);
            assert!(idx % len == 0);
            let slice = &self.fwd_targets[len][idx..idx + len];
            let slice_idx = slice .binary_search_by_key(&(link.lcp() as u32), |x| x.0) .map_or_else(|x| x, |x| x);
            if slice_idx == slice.len() {
                None?;
            }
            let (lcp, target) = slice[slice_idx];
            Link::new(link.source(), link.c(), lcp as usize, target as usize)
        };
        let l = match (k1, k2) {
            (Some(k1), Some(k2)) => Some(std::cmp::min(k1, k2)),
            (Some(k1), None) => Some(k1),
            (None, Some(k2)) => Some(k2),
            (None, None) => None,
        }?;
        if l.source_c() == link.source_c() {
            Some(l)
        } else {
            None
        }
    }
    /// Return the last fwd link <= `link`, but only if it has matching source and character.
    /// FIXME: Currently we implicitly have < instead of <= because the target component is always 0 anyway.
    pub fn fwd_pred(&self, link: Link) -> Option<Link> {
        // eprintln!("PRED OF {link:?}");
        let x = link.key();
        let k1 = self.fwd_ef.pred(x).map(|x| Link::from_key(x.1));
        let k2 = try {
            let &(len, idx) = self.fwd_idx .get(&link.source_c())?;
            let len = len as usize;
            let idx = idx as usize;
            assert!(len > 0);
            assert!(idx % len == 0);
            let slice = &self.fwd_targets[len][idx..idx + len];
            let slice_idx = slice .binary_search_by_key(&(link.lcp() as u32), |x| x.0+1) .map_or_else(|x| x, |x| x+1);
            // eprintln!("slide idx {slice_idx} for link {link:?}");
            if slice_idx == 0 {
                None?;
            }
            let (lcp, target) = slice[slice_idx-1];
            Link::new(link.source(), link.c(), lcp as usize, target as usize)
        };

        let l = std::cmp::max(k1, k2)?;
        if l.source_c() == link.source_c() {
            Some(l)
        } else {
            None
        }
    }
    /// Return an iterator over suffix links >= link.
    pub fn suf_iter_from(&self, link: SuffixLink) -> impl Iterator<Item = SuffixLink> {
        let x = link.key();
        if x <= self.suf_ef.upper_bound() {
            if let Some((_idx, it)) = self.suf_ef.iter_from_succ(x) {
                return Either::Left(
                    it.chain(self.suf_set.iter().copied())
                        .map(SuffixLink::from_key),
                );
            }
        }
        let idx = self.suf_set.binary_search(&x).map_or_else(|x| x, |x| x);
        Either::Right(
            self.suf_set[idx..]
                .iter()
                .copied()
                .map(SuffixLink::from_key),
        )
    }
    // (fwd_ef + suf_ef, fwd_set, suf_set)
    pub fn sizes_gb(&self) -> (f32, f32, f32) {
        (
            gbs(&self.fwd_ef) as f32 + gbs(&self.suf_ef) as f32,
            gbs(&self.fwd_idx) as f32
                + gbs(&self.fwd_targets) as f32
                + gbs(&self.free_targets) as f32,
            gbs(&self.suf_set) as f32,
        )
    }
    pub fn finish(mut self) -> (LinkEf, LinkEf) {
        eprintln!("Merging fwd..");
        let fwd_ef = self.merge_fwd();
        eprintln!("Merging suf..");
        let suf_ef = self.merge_suf();
        (fwd_ef, suf_ef)
    }
    fn merge_fwd(&mut self) -> LinkEf {
        let mut values: Vec<_> = self.fwd_ef.into_iter().collect();
        for (source_c, (len, idx)) in std::mem::take(&mut self.fwd_idx) {
            let len = len as usize;
            let idx = idx as usize;
            let slice = &self.fwd_targets[len][idx..idx + len];
            for &(lcp, target) in slice {
                let (source, c) = Link::unpack_source_c(source_c);
                let link = Link::new(source, c, lcp as usize, target as usize);
                // eprintln!("Link {:?} from idx {idx} len {len}", link);
                values.push(link.key());
            }
        }
        values.voracious_mt_sort(6);
        build_ef(values)
    }

    fn merge_suf(&mut self) -> LinkEf {
        let ef = &self.suf_ef;
        let set = std::mem::take(&mut self.suf_set);
        let n = ef.len() + set.len() - 1;
        let last = ef.upper_bound().max(set.last().copied().unwrap_or(0));
        let mut builder = EliasFanoBuilder::new(n, last);
        // Skip the last EF entry, as it is superseeded by the first list entry.
        for k in ef.into_iter().take(ef.len() - 1).chain(set.into_iter()) {
            builder.push(k);
        }
        builder.build_with_dict()
    }
}
