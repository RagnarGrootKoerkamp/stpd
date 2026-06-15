use itertools::Itertools;
use link::{Link, LinkEf};
use mem_dbg::{MemSize, SizeFlags};
use rayon::prelude::*;
use std::{
    hash::{Hash, Hasher},
    marker::ConstParamTy,
    ops::Range,
    sync::atomic::{AtomicUsize, Ordering},
};
#[allow(unused)]
use sux::{
    dict::EliasFano,
    traits::{Pred, Succ},
};
use voracious_radix_sort::RadixSort;

use crate::{bwt, gbs, lcp::Lcp, sa_and_lcp_cached, T};

mod link;
pub mod storage;
mod stpd;

/// TODO: LexMin, LexMax, CoLexMin, CoLexMax.
#[derive(PartialEq, Eq, ConstParamTy, Debug, Hash)]
pub enum Pi {
    LeftMost,
    RightMost,
}

pub struct JumpIndex<'t, const PI: Pi> {
    pub t: &'t T,
    pub root_anchor: usize,
    pub ef_links: link::LinkEf,

    pub cdawg_nodes: usize,
    pub cdawg_edges: usize,
}

pub struct JumpIndexStats {
    pub num_sampled: usize,
    pub num_sources: usize,
    pub num_source_chars: usize,
    pub num_links: usize,
    pub cdawg_nodes: usize,
    pub cdawg_edges: usize,
}

type ConstructionOutput = (usize, usize, usize, LinkEf);

impl<'t, const PI: Pi> JumpIndex<'t, PI> {
    /// Take an already-built SA, BWT, and LCP.
    ///
    /// These are AsRef so we can give owned objects and drop them as soon as
    /// they are not needed anymore.
    pub fn new(t: &'t T) -> Self {
        let (cdawg_nodes, cdawg_edges, root_anchor, ef_links) = Self::compute_links_cached(t);

        if ef_links.len() < 100 {
            for l in ef_links.iter() {
                eprintln!("link {:?}", Link::from_key(l));
            }
        }

        eprintln!("---");
        eprintln!("final EF size: {}", print_ef(&ef_links));
        eprintln!("---");
        if false {
            {
                eprintln!("splitting.. (drop 1 LCP per (source, c))");
                let (ef_compact, ef_lcp) = link::links_to_compact_ef(&ef_links);
                eprintln!("compact EF size: {}", print_ef(&ef_compact));
                eprintln!("LCP EF size:     {}", print_ef(&ef_lcp));
            }
            eprintln!("---");
            {
                let compact_links = link::compactify(&ef_links);
                eprintln!("compact EF without LCP: {}", print_ef(&compact_links));
            }
            eprintln!("---");
            #[cfg(feature = "mphf")]
            {
                eprintln!("MphfStore.. (dropping (source,c) completely)");
                let store = storage::MphfStore::new(&ef_links);
                eprintln!("{}", store.size());
                eprintln!("---");
            }
        }

        eprintln!(
            "Max LCP: {}",
            ef_links
                .iter()
                .map(|l| Link::from_key(l).lcp())
                .max()
                .unwrap_or(0)
        );

        eprintln!(
            "Average LCP: {:.2}",
            ef_links
                .iter()
                .map(|l| Link::from_key(l).lcp())
                .sum::<usize>() as f32
                / ef_links.len() as f32
        );

        {
            let mut c1 = 0;
            let mut c2 = 0;
            let mut c3 = 0;
            let mut c4 = 0;
            let chunks = ef_links.iter().chunk_by(|l| {
                let l = link::Link::from_key(*l);
                (l.source(), l.c())
            });
            for group in chunks.into_iter().map(|(k, g)| (k, g.count())) {
                c1 += (group.1 > 1) as usize;
                c2 += (group.1 > 2) as usize;
                c3 += (group.1 > 3) as usize;
                c4 += (group.1 > 4) as usize;
            }
            eprintln!("Number of (pos,c) with >1 link: {c1}",);
            eprintln!("Number of (pos,c) with >2 link: {c2}",);
            eprintln!("Number of (pos,c) with >3 link: {c3}",);
            eprintln!("Number of (pos,c) with >4 link: {c4}",);
        }

        JumpIndex {
            t,
            ef_links,
            root_anchor,
            cdawg_nodes,
            cdawg_edges,
        }
    }

    pub fn compute_links_cached<'t2>(t: &'t2 T) -> ConstructionOutput {
        use std::collections::hash_map::DefaultHasher;
        use std::fs::{self, File};
        use std::io::{BufReader, BufWriter};
        use std::path::Path;

        // Hash the input text
        let mut hasher = DefaultHasher::new();
        t.hash(&mut hasher);
        PI.hash(&mut hasher);
        let hash = hasher.finish();

        // Create cache directory if it doesn't exist
        let cache_dir = Path::new("_cache");
        if !cache_dir.exists() {
            fs::create_dir_all(cache_dir).expect("Failed to create _cache directory");
        }

        let cache_file = cache_dir.join(format!("{:x}.ji", hash));

        // Try to load from cache
        if cache_file.exists() {
            eprintln!("Loading from cache: {:?}", cache_file);
            let file = File::open(&cache_file).unwrap();
            let reader = BufReader::new(file);
            let output = bincode::deserialize_from::<_, ConstructionOutput>(reader).unwrap();
            return output;
        }

        // Compute SA and LCP
        let result = Self::compute_links(t);

        // Write to cache
        if t.len() > 100_000_000 {
            eprintln!("Writing to cache: {:?}", cache_file);
            let file = File::create(&cache_file).unwrap();
            let writer = BufWriter::new(file);
            bincode::serialize_into(writer, &result).unwrap();
        }

        result
    }

    /// Returns:
    /// - #cdawg nodes
    /// - #cdawg edges
    /// - root anchor
    /// - EF-encoded links
    fn compute_links<'t2>(t: &'t2 Vec<u8>) -> ConstructionOutput {
        let n = t.len();
        let (sa, lcp) = sa_and_lcp_cached(t);
        let bwt = bwt(t, &sa);
        let link = |anchors: &[usize],
                    lcp: u32,
                    single_run: bool,
                    links: &mut Vec<link::Link>,
                    cdawg_nodes: &mut usize,
                    cdawg_edges: &mut usize|
         -> usize {
            let best = match PI {
                Pi::LeftMost => *anchors.iter().min_by_key(|a| sa[**a]).unwrap(),
                Pi::RightMost => *anchors.iter().max_by_key(|a| sa[**a]).unwrap(),
            };
            // eprintln!("single run: {single_run:?}");
            if single_run || anchors.len() == 1 {
                return best;
            }
            *cdawg_nodes += 1;
            *cdawg_edges += anchors.len();
            for &a in anchors {
                if a == best {
                    continue;
                }
                let text_idx = sa[a] as usize;
                let target = text_idx + lcp as usize;
                // TODO: Why do we need this if statement?
                if target < t.len() {
                    // sa[best] is HOT.
                    let source = sa[best] as usize + lcp as usize;
                    let c = t[target];
                    links.push(link::Link::new(
                        source,
                        c,
                        // co_lcp is HOT.
                        lcp as usize
                            + lcs(&t[..source - lcp as usize], &t[..target - lcp as usize]),
                        target,
                    ));
                    // eprintln!("Link: {:?}", links.last().unwrap());
                }
            }
            best
        };
        let dfs2 = |interval: Range<usize>,
                    links: &mut Vec<link::Link>,
                    cdawg_nodes: &mut usize,
                    cdawg_edges: &mut usize|
         -> usize {
            // eprintln!("DFS2 {interval:?}");
            let mut stack: Vec<(usize, u32, Vec<usize>)> = vec![];

            let mut run_start = interval.start;
            for i in interval.clone() {
                // New run?
                if i > 0 && bwt[i] != bwt[i - 1] {
                    run_start = i;
                }
                // Same run, but one wraps around the text?
                if sa[i] == 0 || (i > 0 && sa[i - 1] == 0) {
                    run_start = i;
                }

                let l = lcp.get(&sa, i);
                // eprintln!("{i} => {l}");

                let mut last_start = i;
                let mut a = i;
                while !stack.is_empty() && (l < stack.last().unwrap().1 || i == interval.end - 1) {
                    let (start, lcp, mut anchors) = stack.pop().unwrap();
                    anchors.push(a);
                    let single_run = run_start <= start;
                    a = link(&anchors, lcp, single_run, links, cdawg_nodes, cdawg_edges);
                    last_start = start;
                    // eprintln!(
                    //     "{start}..={i}: {lcp}  min@{a} {}",
                    //     if single_run { "single run!" } else { "" }
                    // );
                }
                if i == interval.end - 1 {
                    assert!(stack.is_empty());
                    // eprintln!("DFS RETURNS {a}");
                    return a;
                }
                if !stack.is_empty() && l == stack.last().unwrap().1 {
                    stack.last_mut().unwrap().2.push(a);
                } else {
                    // eprintln!("new {last_start}..: {l} first anchor {a}");
                    stack.push((last_start, l, vec![a]));
                }
            }
            unreachable!();

            // assert!(stack.is_empty());
            // links.sort_by_key(Link::key);
            // eprintln!("LINKS");
            // for link in links {
            //     eprintln!("{link:?}");
            // }
        };
        const PREFIX_LCP: u32 = 3;
        let intervals: Vec<usize> = (0..=n)
            .into_par_iter()
            .filter(|&i| i == 0 || lcp.get(&sa, i - 1) <= PREFIX_LCP)
            .collect();
        eprintln!("Run on {} intervals!", intervals.len());
        let done = std::sync::atomic::AtomicUsize::new(0);
        let total = std::sync::atomic::AtomicUsize::new(0);
        let ef_total = std::sync::atomic::AtomicUsize::new(0);
        let cdawg_nodes = AtomicUsize::new(0);
        let cdawg_edges = AtomicUsize::new(0);
        let (anchors, mut dfs_results): (Vec<usize>, Vec<link::BareEf>) = intervals
            .par_array_windows()
            // .array_windows()
            .map(|&[start, end]| {
                // eprintln!("sa[{start}..{end}]");
                let mut links = vec![];
                let mut cn = 0;
                let mut ce = 0;
                let a = dfs2(start..end, &mut links, &mut cn, &mut ce);
                let links_ef = link::links_to_ef(links);

                let _done = done.fetch_add(1, Ordering::Relaxed) + 1;
                let _total = total.fetch_add(links_ef.len(), Ordering::Relaxed) + links_ef.len();
                let ef_size = MemSize::mem_size(&links_ef, SizeFlags::default());
                let _ef_total = ef_total.fetch_add(ef_size, Ordering::Relaxed) + ef_size;
                cdawg_nodes.fetch_add(cn, Ordering::Relaxed);
                cdawg_edges.fetch_add(ce, Ordering::Relaxed);

                // eprintln!(
                //     "{done:>2}: Collected {} links, {total} total, EF {:.3} GB (total EF {:.3} GB; total flat {:.3})",
                //     links_ef.len(),
                //     ef_size as f32 / 1e9,
                //     ef_total as f32 / 1e9,
                //     (total * std::mem::size_of::<Link>()) as f32 / 1e9,
                // );

                (a, links_ef)
            })
            .unzip();

        let mut cdawg_nodes = cdawg_nodes.into_inner();
        let mut cdawg_edges = cdawg_edges.into_inner();

        let mut links: Vec<link::Link> = vec![];
        let root_anchor;
        {
            let mut intervals = intervals;
            let mut anchors = anchors;
            // eprintln!("Intervals: {intervals:?}");
            // eprintln!("Anchors:   {anchors:?}");
            for cur_lcp in (1..=PREFIX_LCP + 0).rev() {
                // eprintln!("cur lcp: {cur_lcp}");
                let mut new_intervals = vec![0];
                let mut new_anchors = vec![];

                let mut start_idx = 0;
                // for start_idx in 0..intervals.len() - 1 {
                for end_idx in 1..intervals.len() {
                    // let end_idx = start_idx + 1;
                    let start = intervals[start_idx];
                    let end = intervals[end_idx];
                    if lcp.get(&sa, end - 1) < cur_lcp {
                        // eprintln!("intervals[{start_idx}..{end_idx}]");
                        // eprintln!("sa[{start}..{end}]");
                        new_intervals.push(end);
                        new_anchors.push(link(
                            &anchors[start_idx..end_idx],
                            cur_lcp,
                            bwt[start..end].iter().all_equal()
                                && sa[start..end].iter().all(|&x| x != 0),
                            &mut links,
                            &mut cdawg_nodes,
                            &mut cdawg_edges,
                        ));
                        start_idx = end_idx;
                    }
                }
                anchors = new_anchors;
                intervals = new_intervals;
                // eprintln!("Intervals: {intervals:?}");
                // eprintln!("Anchors:   {anchors:?}");
            }
            let idx_of_min = link(
                &anchors,
                0,
                bwt.iter().all_equal(),
                &mut links,
                &mut cdawg_nodes,
                &mut cdawg_edges,
            );
            eprintln!("Idx of min: {idx_of_min}");
            root_anchor = sa[idx_of_min] as usize;
            eprintln!("Root anchor: {root_anchor}");
            match PI {
                Pi::LeftMost => {
                    assert_eq!(root_anchor, 0);
                }
                Pi::RightMost => {
                    assert_eq!(root_anchor, sa.len() - 1);
                }
            }
        }
        {
            // Collect links from top layers.
            let mut links = links.into_iter().map(|l| l.key()).collect_vec();
            links.sort();
            links.dedup();
            dfs_results.push(link::BareEf::from(links));
        }
        drop(bwt);
        drop(lcp);
        drop(sa);
        let ef_links = Self::merge_links(dfs_results);
        (cdawg_nodes, cdawg_edges, root_anchor, ef_links)
    }

    fn merge_links(dfs_results: Vec<EliasFano<u128>>) -> link::LinkEf {
        let num_vals = dfs_results.iter().map(|x| x.len()).sum();

        eprintln!(
            "Total EF size before merging: {:.3} GB",
            dfs_results
                .iter()
                .map(|ef| mem_dbg::MemSize::mem_size(ef, mem_dbg::SizeFlags::default()))
                .sum::<usize>() as f32
                / 1e9
        );

        eprintln!("Select pivots");
        const NUM_PIVOTS: usize = 80;
        const OVERSAMPLING: usize = 100;
        // Oversample 100x to smoothen the distribution.
        let mut pivot_idxs = (0..NUM_PIVOTS * OVERSAMPLING)
            .map(|_| rand::random_range(0..num_vals))
            .collect_vec();
        pivot_idxs.sort_unstable();
        pivot_idxs.dedup();

        let mut pivots = vec![];
        {
            let mut i = 0;
            let mut pivot_idx_iter = pivot_idxs.iter();
            let mut next_pivot = *pivot_idx_iter.next().unwrap();
            'outer: for ef in &dfs_results {
                for x in ef.iter() {
                    if i == next_pivot {
                        pivots.push(x);
                        next_pivot = match pivot_idx_iter.next() {
                            Some(&idx) => idx,
                            None => break 'outer,
                        };
                    }
                    i += 1;
                }
            }
            pivots.sort_unstable();
            // Keep only every 100'th pivot
            pivots = pivots
                .into_iter()
                .enumerate()
                .filter(|(i, _)| i % OVERSAMPLING == 0)
                .map(|(_, x)| x)
                .collect_vec();
            pivots.insert(0, 0);
            pivots.push(u128::MAX);
            pivots.dedup();
        }

        let efs_per_pivot: Vec<std::sync::Mutex<Vec<link::BareEf>>> = (0..pivots.len())
            .map(|_| std::sync::Mutex::new(vec![]))
            .collect();
        let max = std::sync::Mutex::new(0u128);
        eprintln!("Partitioning EFs");
        dfs_results.into_par_iter().for_each(|ef| {
            if ef.len() == 0 {
                return;
            }
            let vals = ef.iter().collect_vec();
            {
                let mut max = max.lock().unwrap();
                *max = (*max).max(*vals.last().unwrap());
            }
            let splits = pivots
                .iter()
                .map(|&p| vals.partition_point(|&x| x < p))
                .collect_vec();
            for i in 0..pivots.len() - 1 {
                let bare_ef = link::BareEf::from(&vals[splits[i]..splits[i + 1]]);
                efs_per_pivot[i].lock().unwrap().push(bare_ef);
            }
        });
        let efs_per_pivot = efs_per_pivot
            .into_iter()
            .map(|m| m.into_inner().unwrap())
            .collect_vec();
        let max = max.into_inner().unwrap();
        eprintln!(
            "Total EF size after partitioning: {:.3} GB",
            efs_per_pivot.iter().flatten().map(crate::gbs).sum::<f32>()
        );

        eprintln!("Build an EF for each part");
        let part_efs: Vec<(u128, link::BareEf)> = efs_per_pivot
            .into_par_iter()
            .enumerate()
            .map(|(_i, efs)| {
                let mut vals = vec![];
                // eprintln!(
                //     "{i}: Merging {} EFs of total len {} total size {:.3}",
                //     efs.len(),
                //     efs.iter().map(|ef| ef.len()).sum::<usize>(),
                //     efs.iter().map(crate::gbs).sum::<f32>()
                // );
                for ef in efs {
                    vals.extend(ef.iter());
                }
                // eprintln!(
                //     "{i}: vals size: {:.3} GB",
                //     std::mem::size_of_val(vals.as_slice()) as f32 / 1e9
                // );
                vals.voracious_sort();
                // vals.sort_unstable();
                vals.dedup();
                // eprintln!(
                //     "{i}: vals size after dedup: {:.3} GB",
                //     std::mem::size_of_val(vals.as_slice()) as f32 / 1e9
                // );
                let min = *vals.first().unwrap_or(&0);
                for x in &mut *vals {
                    *x -= min;
                }
                let out = link::BareEf::from(vals);
                // eprintln!(
                //     "{i}: output EF size: {:.3} GB",
                //     mem_dbg::MemSize::mem_size(&out, mem_dbg::SizeFlags::default()) as f32
                //         / 1e9
                // );
                (min, out)
            })
            .collect();
        eprintln!(
            "Total EF size after building parts: {:.3} GB",
            part_efs
                .iter()
                .map(|(_min, ef)| mem_dbg::MemSize::mem_size(ef, mem_dbg::SizeFlags::default()))
                .sum::<usize>() as f32
                / 1e9
        );

        let n = part_efs.iter().map(|ef| ef.1.len()).sum::<usize>();
        // eprintln!("Merge part EFs. Dedupped to {n} links");
        let mut ef_builder = sux::dict::elias_fano::EliasFanoBuilder::<u128>::new(n, max);
        for (min, part_ef) in part_efs {
            // eprintln!(
            //     "Extending by {} elems size {}",
            //     part_ef.len(),
            //     mem_dbg::MemSize::mem_size(&part_ef, mem_dbg::SizeFlags::default()) as f32
            //         / 1e9
            // );
            ef_builder.extend(part_ef.iter().map(|x| x + min));
        }
        ef_builder.build_with_dict()
    }

    pub fn stats(&self) -> JumpIndexStats {
        let mut stpd_samples: Vec<usize> = self
            .ef_links
            .iter()
            .map(|l| link::Link::from_key(l).target())
            .collect();
        stpd_samples.voracious_mt_sort(12);
        stpd_samples.dedup();
        JumpIndexStats {
            num_sampled: stpd_samples.len(),
            num_sources: 1 + self
                .ef_links
                .iter()
                .tuple_windows()
                .filter(|&(a, b)| {
                    link::Link::from_key(a).source() != link::Link::from_key(b).source()
                })
                .count(),
            num_source_chars: 1 + self
                .ef_links
                .iter()
                .tuple_windows()
                .filter(|&(a, b)| {
                    link::Link::from_key(a).source_c() != link::Link::from_key(b).source_c()
                })
                .count(),
            num_links: 1 + self
                .ef_links
                .iter()
                .tuple_windows()
                .filter(|(a, b)| a != b)
                .count(),
            cdawg_nodes: self.cdawg_nodes,
            cdawg_edges: self.cdawg_edges,
        }
    }

    pub fn inspect_links(&self) {
        eprintln!("Links: {}", self.ef_links.len());
        if self.ef_links.len() < 10000 {
            // for i in 0..self.ef_links.len() {
            //     eprintln!("{i:>8} {:?}", self.ef_links[i]);
            // }
        } else {
            // for i in (0..self.ef_links.len()).step_by(100000) {
            //     for i in i..i + 100 {
            //         eprintln!("{i:>8} {:?}", self.ef_links[i]);
            //     }
            //     eprintln!("---");
            // }
        }
        // eprintln!("---");

        // eprintln!("Samples: {}", self.stpd_samples.len());
        // if self.stpd_samples.len() < 30000 {
        //     for i in 0..self.stpd_samples.len() {
        //         eprintln!("{i:>8} {:>7}", self.stpd_samples[i]);
        //     }
        // } else {
        //     for i in (0..self.stpd_samples.len()).step_by(100000) {
        //         for i in i..i + 100 {
        //             eprintln!("{i:>8} {:>7}", self.stpd_samples[i]);
        //         }
        //         eprintln!("---");
        //     }
        // }
    }

    pub fn space(&self) {
        eprintln!(
            "jump index   {:.3} GB",
            mem_dbg::MemSize::mem_size(&self.ef_links, mem_dbg::SizeFlags::default()) as f32 / 1e9
        );
    }

    /// Returns the position in the text of the longest prefix of the pattern that matches.
    /// Also returns the total number of jumps.
    pub fn map_jump(&self, pattern: &[u8]) -> (Range<usize>, usize) {
        let mut pos = self.root_anchor;
        let mut jumps = 0;
        for (i, &c) in pattern.iter().enumerate() {
            if pos < self.t.len() && self.t[pos] == c {
                pos += 1;
                continue;
            }
            // if pos < self.t.len() {
            //     eprintln!("{i}: mismatch at {pos}: got {} wanted {c}", self.t[pos]);
            // } else {
            //     eprintln!("{i}: mismatch at {pos}: got end of text wanted {c}",);
            // }

            // find the first link at (pos, c) with LCP >= i.
            let Some((_idx, link)) = self.ef_links.succ(&link::Link::new(pos, c, i, 0).key())
            else {
                // key was larger than last link
                return (pos - i..pos, jumps);
            };
            let link = link::Link::from_key(link);
            if link.source() == pos && link.c() as u8 == c {
                // eprintln!("link: {link:?}");
                jumps += 1;
                pos = link.target() + 1;
            } else {
                // eprintln!("no link found; next is {link:?}");
                // if let Some((_idx, link)) = self.ef_links.pred(&link::Link::new(pos, c, i, 0).key())
                // {
                //     let link = link::Link::from_key(link);
                //     eprintln!("prev link: {link:?}");
                // }
                return (pos - i..pos, jumps);
            }
        }
        // eprintln!("Pattern found at {}", pos-pattern.len());
        (pos - pattern.len()..pos, jumps)
    }

    /// Relative Lempel-Ziv.
    /// Repeatedly greedily matches a longest prefix of the remaining pattern.
    pub fn map_rlz(&self, pattern: &[u8]) -> (usize, usize) {
        let mut start = 0;
        let mut parts = 0;
        let mut jumps = 0;
        while start < pattern.len() {
            let (range, js) = self.map_jump(&pattern[start..]);
            start += range.len();
            parts += 1;
            jumps += js;
        }
        assert_eq!(start, pattern.len());
        (parts, jumps)
    }

    /// Take a bunch of random substrings and map them against the text.
    pub fn test_map(&self) {
        let cnt = 10000000;
        let len = 1..self.t.len().min(5000);
        for _ in 0..cnt {
            let len = rand::random_range(len.clone());
            let i = rand::random_range(0..=self.t.len() - len);
            let j = i + len;

            let pattern = &self.t[i..j];
            // eprintln!("Searching pattern T[{i}..{i}+{len}]");
            let p1 = self.map_jump(pattern).0;
            assert!(p1.len() == pattern.len(), "substring {i}..{j} not found");
            let pos = p1.start;
            match PI {
                Pi::LeftMost => assert!(pos <= i, "substring {i}..{j} found at pos {pos}"),
                Pi::RightMost => assert!(pos >= i, "substring {i}..{j} found at pos {pos}"),
            }
            // eprintln!("substring {i}..{j} found at pos {pos}");
        }
    }

    /// Take a bunch of random substrings and map them against the text.
    pub fn bench_rlz(&self) {
        let cnt = 100000;
        let len = 1..5000.min(self.t.len());

        for rate in [0.01, 0.001] {
            let mut patterns = vec![];
            for _it in 0..cnt {
                let len = rand::random_range(len.clone());
                let i = rand::random_range(0..=self.t.len() - len);
                let j = i + len;
                let mut pattern = self.t[i..j].to_vec();
                // randomly permute 1% of values.
                for _ in 0..(len as f32 * rate) as usize {
                    let idx = rand::random_range(0..len);
                    pattern[idx] = (pattern[idx] + rand::random::<u8>() % 3) % 4;
                }
                patterns.push(pattern);
            }
            let start = std::time::Instant::now();
            let mut parts = 0;
            let mut jumps = 0;
            for (_it, pattern) in patterns.iter().enumerate() {
                let (ps, js) = self.map_rlz(&pattern);
                parts += ps;
                jumps += js;
            }
            let dur = start.elapsed();
            let avg_part_len =
                patterns.iter().map(|p| p.len()).sum::<usize>() as f32 / parts as f32;
            let avg_jump_dist =
                patterns.iter().map(|p| p.len()).sum::<usize>() as f32 / jumps as f32;
            let avg_jumps_per_part = jumps as f32 / parts as f32;
            eprintln!(
                "RLZ: {:.3} for {} reads of average length 2500 with {:4.2}% errors. Avg {:6.0} reads/sec. Avg part len: {avg_part_len:6.1}, avg jump dist: {avg_jump_dist:6.1}, avg jumps/part: {avg_jumps_per_part:6.2}",
            dur.as_secs_f32(),
            cnt,
                rate*100.,
            cnt as f32 / dur.as_secs_f32(),
        );
        }
    }
}

fn lcs(a: &[u8], b: &[u8]) -> usize {
    let min = a.len().min(b.len());
    for i in 0..min {
        if a[a.len() - 1 - i] != b[b.len() - 1 - i] {
            return i;
        }
    }
    min
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::print;
    use crate::strings::relative;

    #[test]
    fn fuzz() {
        test_direction::<{ Pi::LeftMost }>();
        test_direction::<{ Pi::RightMost }>();
    }
    fn test_direction<const PI: Pi>() {
        eprintln!("--- TESTING {PI:?} ---");
        // let mut t1 = std::time::Duration::default();
        let mut t2 = std::time::Duration::default();
        for (len, repeats, r) in [
            (10, 1, 0.1),
            (10, 2, 0.1),
            (10, 3, 0.1),
            (10, 4, 0.1),
            (10, 10, 0.1),
            (100, 1, 0.1),
            (100, 10, 0.1),
            (100, 100, 0.05),
            (100, 1000, 0.05),
            (1000, 100, 0.05),
            (10000, 10, 0.05),
        ] {
            for _ in 0..10 {
                let t = relative(len, 4, repeats, r).1;
                eprintln!("text: {}", print(&t));

                eprintln!("building for {len}x{repeats} at {r}..");
                let ji = JumpIndex::<PI>::new(&t);

                let maxlen = t.len().min(1000);
                eprintln!("querying..");
                for _id in 0..100000 {
                    // FIXME: Include the empty suffix in the suffix array.
                    let len = rand::random_range(1..=maxlen);

                    let pos = rand::random_range(0..=t.len() - len);
                    let pattern = &t[pos..pos + len];

                    eprintln!(
                        "pattern for {PI:?}: T[{pos}..{pos}+{len}] = {}",
                        print(pattern)
                    );

                    // let s = std::time::Instant::now();
                    // let p1 = ji.map_stpd(pattern);
                    // eprintln!("p1: {p1:?}");
                    // t1 += s.elapsed();
                    let s = std::time::Instant::now();
                    let p2 = ji.map_jump(pattern).0;
                    t2 += s.elapsed();
                    eprintln!("p2: {p2:?}");
                    assert_eq!(p2.len(), len, "Did not match the full pattern!");

                    match PI {
                        Pi::LeftMost => {
                            assert!(
                                p2.start <= pos,
                                "substring {pos}..{pos}+{len} found at pos {p2:?}"
                            )
                        }
                        Pi::RightMost => {
                            assert!(
                                p2.start >= pos,
                                "substring {pos}..{pos}+{len} found at pos {p2:?}"
                            )
                        }
                    }

                    // let p1 = p1.unwrap();
                    // let p2 = p2;

                    // assert_eq!(&t[p1..p1 + len], pattern);
                    if &t[p2.clone()] != pattern {
                        eprintln!("Pattern T[{pos}..{pos}+{len}] does not match text T[{p2:?}] for JI<{PI:?}>!");
                        eprintln!("pattern: {}", print(pattern));
                        eprintln!("text:    {}", print(&t[p2]));
                        panic!();
                    }
                    // assert_eq!(p1, p2);
                }
            }
            // eprintln!("STPD: {t1:?}");
            eprintln!("JI:   {t2:?}");
        }
    }

    fn test_one(t: &Vec<u8>, pos: usize, len: usize) {
        const PI: Pi = Pi::LeftMost;
        eprintln!("text: {}", print(&t));

        let ji = JumpIndex::<PI>::new(t);
        let pattern = &t[pos..pos + len];

        eprintln!(
            "pattern for {PI:?}: T[{pos}..{pos}+{len}] = {}",
            print(pattern)
        );

        let p2 = ji.map_jump(pattern).0;
        assert_eq!(p2.len(), len, "Did not match the full pattern!");

        if &t[p2.clone()] != pattern {
            eprintln!(
                "Pattern T[{pos}..{pos}+{len}] does not match text T[{p2:?}] for JI<{PI:?}>!"
            );
            eprintln!("pattern: {}", print(pattern));
            eprintln!("text:    {}", print(&t[p2]));
            panic!();
        }
    }

    /// Error: Missing link for CCA at pos 2 because CCA and CCB are a BTW-run
    /// but the character preceding CCB wraps around the text.
    #[test]
    fn failure_one() {
        test_one(&b"CCBDDACADDCCADDACBDD".to_vec(), 10, 9);
    }

    /// Like failure_one, but the BWT-run for the first character is long and in the DFS2 recursion.
    #[test]
    fn failure_two() {
        test_one(&b"DCACCCABDCDCACCCABDCDCCCCCABAC".to_vec(), 19, 7);
    }

    /// Like failure_one and failure_two, to catch that both (x, 0) and (0, x)
    /// are run-breaking in the suffix array.
    #[test]
    fn failure_three() {
        test_one(&b"CCBDDACADDCCADDACBDDCACADDCCADDACBDD".to_vec(), 16, 12);
    }
}

fn print_ef<V, H, L>(ef: &EliasFano<V, H, L>) -> String
where
    EliasFano<V, H, L>: MemSize,
{
    let n = ef.len();
    let l = ef.num_lower_bits();
    format!("{:.3} GB : {n} * (2 + {l}) bits", gbs(ef))
}

// fn build_ef(vals: &[u128]) -> EliasFanoConcurrentBuilder<u64> {
//     let builder = EliasFanoConcurrentBuilder::new(vals.len(), *vals.last().unwrap());
//     vals.par_iter().enumerate().for_each(|(i, &v)| {
//         builder.set(i, v);
//     });
//     builder
// }
