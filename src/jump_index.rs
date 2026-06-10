use itertools::Itertools;
use std::{
    cmp::Ordering::{Greater, Less},
    marker::Sync,
};
use sux::{bits::BitVec, dict::EliasFano, rank_sel::SelectZeroAdaptConst, traits::Succ};
use voracious_radix_sort::RadixSort;

use crate::{
    bwt,
    lcp::CompactLcp,
    longest_common_prefix,
    rmq::{self, Rmq},
    sa_and_lcp,
    stpd::cmp_colex,
    SaElem, SA, T,
};

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct Link {
    data: u128,
    // source: usize,
    // c: u8,
    // lcp: usize,
    // target: usize,
}

const SOURCE_BITS: u32 = 31;
const C_BITS: u32 = 8; // TODO: Reduce to 2
const LCP_BITS: u32 = 22; // enough for reference genome
const TARGET_BITS: u32 = 31; // enough for reference genome
const LINK_BITS: u32 = SOURCE_BITS + C_BITS + LCP_BITS + TARGET_BITS;

type LinkEf = EliasFano<u128, SelectZeroAdaptConst<BitVec<Box<[usize]>>, Box<[usize]>, 12, 3>>;
type BareEf = EliasFano<u128, BitVec<Box<[usize]>>>;

impl Link {
    const MAX: u128 = (1 << LINK_BITS) - 1;
    fn from_key(data: u128) -> Self {
        Self { data }
    }
    fn new(source: usize, c: u8, lcp: usize, target: usize) -> Self {
        assert!(LINK_BITS <= 128);
        assert!(
            source < (1 << SOURCE_BITS) as usize,
            "link {source},{c} -> {lcp},{target}"
        );
        assert!(
            (c as usize) < (1 << C_BITS),
            "link {source},{c} -> {lcp},{target}"
        );
        assert!(
            lcp < (1 << LCP_BITS) as usize,
            "link {source},{c} -> {lcp},{target}"
        );
        assert!(
            target < (1 << TARGET_BITS) as usize,
            "link {source},{c} -> {lcp},{target}"
        );
        let data = ((source as u128) << (C_BITS + LCP_BITS + TARGET_BITS))
            | ((c as u128) << (LCP_BITS + TARGET_BITS))
            | ((lcp as u128) << TARGET_BITS)
            | (target as u128);
        Self { data }
    }
    fn key(&self) -> u128 {
        self.data
    }
    fn source(&self) -> usize {
        ((self.key() >> (C_BITS + LCP_BITS + TARGET_BITS)) & ((1 << SOURCE_BITS) - 1)) as usize
    }
    fn source_c(&self) -> usize {
        ((self.key() >> (LCP_BITS + TARGET_BITS)) & ((1 << (SOURCE_BITS + C_BITS)) - 1)) as usize
    }
    fn c(&self) -> u8 {
        ((self.key() >> (LCP_BITS + TARGET_BITS)) & ((1 << C_BITS) - 1)) as u8
    }
    fn lcp(&self) -> usize {
        ((self.key() >> TARGET_BITS) & ((1 << LCP_BITS) - 1)) as usize
    }
    fn target(&self) -> usize {
        (self.key() & ((1 << TARGET_BITS) - 1)) as usize
    }

    /// Store all values 'mirrored': `u-x` for x from large to small,
    /// where `u` is the largest element.
    fn links_to_ef(links: Vec<Link>) -> BareEf {
        let mut links = links.into_iter().map(|l| l.key()).collect_vec();
        links.voracious_sort();
        links.dedup();
        BareEf::from(links)
    }
}

pub struct JumpIndex<TR: AsRef<T>> {
    pub t: TR,
    pub stpd_samples: Vec<usize>,
    pub stpd_pi: Vec<u64>,
    pub stpd_rmq: rmq::BlockRmq<u64, 64>,
    // TODO: Predecessor structure
    pub ef_links: LinkEf,
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

impl<TR: AsRef<T> + Sync> JumpIndex<TR> {
    pub fn new(t: TR) -> Self {
        let (sa, lcp) = sa_and_lcp(t.as_ref());
        let bwt = &bwt(t.as_ref(), &sa);
        // let pi = (0..t.as_ref().len()).collect_vec();
        Self::new2(t, sa, bwt, &lcp, &vec![])
    }
}
impl<TR: AsRef<T> + Sync> JumpIndex<TR> {
    pub fn new2<SAR: AsRef<SA> + Sync>(t: TR, sa: SAR, bwt: &T, lcp: &CompactLcp, pi: &SA) -> Self {
        struct State<'a, T, SA> {
            t: T,
            bwt: &'a Vec<u8>,
            sa: &'a SA,
            lcp: &'a CompactLcp,
            // run_boundaries: Vec<SaElem>,
            lcp_rmq: rmq::BlockRmq<crate::LcpElem, 32>,
            permuted_pi: &'a SA,
            pi_rmq: rmq::BlockRmq<SaElem, 128>,
        }

        impl<'a, TR: AsRef<T>, SAR: AsRef<SA>> State<'a, TR, SAR> {
            fn split(
                &self,
                interval: std::ops::Range<usize>,
            ) -> Option<(bool, usize, usize, Vec<std::ops::Range<usize>>)> {
                if interval.len() <= 1 {
                    return None;
                }
                // let idx = self
                //     .run_boundaries
                //     .binary_search(&(interval.start as SaElem))
                //     .unwrap_or_else(|x| x);
                // let single_run = idx == self.run_boundaries.len()
                //     || self.run_boundaries[idx] >= interval.end as SaElem - 1;
                // FIXME: Is this too slow? Use EF on run boundaries instead?
                let single_run = self.bwt[interval.clone()]
                    .iter()
                    .all(|&c| c == self.bwt[interval.start]);
                if single_run {
                    return None;
                }

                let anchor_idx = self
                    .pi_rmq
                    .query(
                        self.permuted_pi.as_ref().as_slice(),
                        interval.start,
                        interval.end - 1,
                    )
                    .1;
                // let _anchor_pos = self.sa.as_ref()[anchor_idx];
                // eprintln!("anchor pos: {anchor_pos}");
                let mut done_intervals = vec![];
                let mut wip_intervals = vec![interval.clone()];
                // FIXME: Keep track of LCP so far in DFS on suffix tree.
                let lcp = longest_common_prefix(
                    &self.t.as_ref()[self.sa.as_ref()[interval.start] as usize..],
                    &self.t.as_ref()[self.sa.as_ref()[interval.end - 1] as usize..],
                );
                // let lcp = self.lcp.get(
                //     self.sa.as_ref(),
                //     self.lcp_rmq
                //         .query(&self.lcp, interval.start, interval.end - 2)
                //         .1,
                // );
                // eprintln!(
                //     "node for {}",
                //     crate::print(&self.t.as_ref()[anchor_pos..anchor_pos + lcp as usize])
                // );

                while let Some(interval) = wip_intervals.try_remove(0) {
                    if interval.len() <= 1 {
                        done_intervals.push(interval);
                        continue;
                    }

                    let split_pos = self
                        .lcp_rmq
                        .query(
                            (self.sa.as_ref(), self.lcp),
                            interval.start,
                            interval.end - 2,
                        )
                        .1
                        + 1;
                    // FIXME: Keep track of LCP so far in DFS on suffix tree.
                    let new_lcp = longest_common_prefix(
                        &self.t.as_ref()[self.sa.as_ref()[split_pos - 1] as usize..],
                        &self.t.as_ref()[self.sa.as_ref()[split_pos] as usize..],
                    );
                    // let new_lcp = self.lcp.get(self.sa.as_ref(), split_pos - 1);
                    if new_lcp > lcp {
                        done_intervals.push(interval);
                        continue;
                    }
                    assert!(new_lcp == lcp);
                    wip_intervals.push(interval.start..split_pos);
                    wip_intervals.push(split_pos..interval.end);
                }
                Some((single_run, anchor_idx, lcp as usize, done_intervals))
            }

            fn node_output(
                &self,
                anchor_pos: usize,
                lcp: usize,
                done_intervals: &[std::ops::Range<usize>],
                _sampled: &mut Vec<usize>,
                links: &mut Vec<Link>,
            ) {
                for x in done_intervals {
                    if !x.contains(&anchor_pos) {
                        let secondary_anchor_pos = self
                            .pi_rmq
                            .query(self.permuted_pi.as_ref().as_slice(), x.start, x.end - 1)
                            .1;
                        let text_idx = self.sa.as_ref()[secondary_anchor_pos] as usize;
                        let target = text_idx + lcp as usize;
                        if target < self.t.as_ref().len() {
                            // NOTE: STPD samples are skipped for now.
                            // sampled.push(target);
                            let source = self.sa.as_ref()[anchor_pos] as usize + lcp as usize;
                            let c = self.t.as_ref()[target];
                            links.push(Link::new(
                                source,
                                c,
                                co_lcp(&self.t.as_ref()[..source], &self.t.as_ref()[..target]),
                                target,
                            ));
                        }
                    }
                }
            }

            /// Returns a mask of BWT chars in interval.
            fn dfs(
                &self,
                interval: std::ops::Range<usize>,
                sampled: &mut Vec<usize>,
                links: &mut Vec<Link>,
                cdawg_nodes: &mut usize,
                cdawg_edges: &mut usize,
            ) {
                assert!(interval.len() > 0);
                let Some((single_run, anchor_pos, lcp, done_intervals)) =
                    self.split(interval.clone())
                else {
                    return;
                };
                assert!(!single_run);
                assert!(done_intervals.len() > 1);

                self.node_output(anchor_pos, lcp, &done_intervals, sampled, links);
                *cdawg_nodes += 1;
                *cdawg_edges += done_intervals.len();

                for x in &done_intervals {
                    self.dfs(x.clone(), sampled, links, cdawg_nodes, cdawg_edges);
                }
            }

            fn collect_work(
                &self,
                interval: std::ops::Range<usize>,
                sampled: &mut Vec<usize>,
                links: &mut Vec<Link>,
                work_queue: &mut Vec<std::ops::Range<usize>>,
                cdawg_nodes: &mut usize,
                cdawg_edges: &mut usize,
            ) {
                const TARGET_CHUNKS: usize = 128;
                if interval.len() < self.t.as_ref().len().div_ceil(TARGET_CHUNKS).max(1_000_000) {
                    work_queue.push(interval);
                    return;
                }
                let Some((single_run, anchor_pos, lcp, done_intervals)) =
                    self.split(interval.clone())
                else {
                    return;
                };
                // eprintln!(
                //     "Split interval of len {} into {} parts",
                //     interval.len(),
                //     done_intervals.len()
                // );
                assert!(!single_run);
                assert!(done_intervals.len() > 1);
                *cdawg_nodes += 1;
                *cdawg_edges += done_intervals.len();

                self.node_output(anchor_pos, lcp, &done_intervals, sampled, links);
                for x in done_intervals {
                    self.collect_work(x, sampled, links, work_queue, cdawg_nodes, cdawg_edges);
                }
            }
        }

        // let run_boundaries: Vec<SaElem> = (0..t.as_ref().len() - 1)
        //     .into_par_iter()
        //     .filter(|&i| bwt[i] != bwt[i + 1])
        //     .map(|i| i as SaElem)
        //     .collect();
        // eprintln!(
        //     "Run boundaries: {:.3} GB",
        //     std::mem::size_of_val(run_boundaries.as_slice()) as f32 / 1e9
        // );
        // let run_boundaries = BTreeSet::from_iter(run_boundaries);
        // let s = mem_dbg::MemSize::mem_size(&run_boundaries, SizeFlags::default());
        // eprintln!("Run boundaries set: {:.3} GB", s as f32 / 1e9);

        // empty indicates pi=identity and permuted_pi = sa.
        let permuted_pi: &SA = if pi.is_empty() {
            &sa.as_ref()
        } else {
            eprintln!("permuting pi..");
            &sa.as_ref().par_iter().map(|&i| pi[i as usize]).collect()
        };

        // eprintln!();
        // eprintln!("sa:  {:?}", sa.as_ref());
        // eprintln!("lcp: {:?}", lcp.as_ref());
        // eprintln!("ppi: {:?}", permuted_pi);
        // eprintln!("run boundaries: {:?}", run_boundaries);

        use rmq::Rmq as _;
        let lcp_rmq = rmq::BlockRmq::build((sa.as_ref(), lcp));
        eprintln!("lcp rmq: {:.3} GB", lcp_rmq.space() as f32 / 1e9);
        let pi_rmq = rmq::BlockRmq::build(permuted_pi.as_slice());
        eprintln!("pi  rmq: {:.3} GB", pi_rmq.space() as f32 / 1e9);
        let state = State {
            t,
            bwt,
            sa: sa.as_ref(),
            // run_boundaries,
            lcp_rmq,
            lcp,
            pi_rmq,
            permuted_pi: &permuted_pi,
        };

        let mut stpd_samples = vec![];
        let mut links = vec![];
        let mut work_queue = vec![];
        let mut cdawg_nodes = 1; // the final node
        let mut cdawg_edges = 0;
        eprintln!("collecting work..");
        state.collect_work(
            0..state.t.as_ref().len(),
            &mut stpd_samples,
            &mut links,
            &mut work_queue,
            &mut cdawg_nodes,
            &mut cdawg_edges,
        );

        eprintln!("Run on {} intervals!", work_queue.len());
        let done = std::sync::atomic::AtomicUsize::new(0);
        let total = std::sync::atomic::AtomicUsize::new(0);
        let ef_total = std::sync::atomic::AtomicUsize::new(0);
        use rayon::prelude::*;
        let mut child_results: Vec<(Vec<usize>, BareEf, usize, usize)> = work_queue
            .into_par_iter()
            .map(|interval| {
                let mut sampled = vec![];
                let mut links = vec![];
                let mut cdawg_nodes = 0;
                let mut cdawg_edges = 0;
                state.dfs(
                    interval,
                    &mut sampled,
                    &mut links,
                    &mut cdawg_nodes,
                    &mut cdawg_edges,
                );
                let links_ef = Link::links_to_ef(links);

                let done = done.fetch_add(1, std::sync::atomic::Ordering::SeqCst) + 1;
                let total =
                    total.fetch_add(links_ef.len(), std::sync::atomic::Ordering::SeqCst) + links_ef.len();
                let ef_size = mem_dbg::MemSize::mem_size(&links_ef, mem_dbg::SizeFlags::default());
                let ef_total =
                    ef_total.fetch_add(ef_size, std::sync::atomic::Ordering::SeqCst) + ef_size;

                eprintln!(
                    "{done:>2}: Collected {} links, {total} total, EF {:.3} GB (total EF {:.3} GB; total flat {:.3})",
                    links_ef.len(),
                    ef_size as f32 / 1e9,
                    ef_total as f32 / 1e9,
                    (total * std::mem::size_of::<Link>()) as f32 / 1e9,
                );

                (sampled, links_ef, cdawg_nodes, cdawg_edges)
            })
            .collect();

        {
            let mut links = links.into_iter().map(|l| l.key()).collect_vec();
            links.sort();
            links.dedup();
            child_results.push((vec![], BareEf::from(links), 0, 0));
        }

        let mut num_vals = 0;
        for (_s, ef, cn, ce) in &child_results {
            num_vals += ef.len();
            // stpd_samples.extend(s);
            cdawg_nodes += cn;
            cdawg_edges += ce;
        }

        let ef_links = {
            eprintln!(
                "Total EF size before merging: {:.3} GB",
                child_results
                    .iter()
                    .map(|(_s, ef, _cn, _ce)| mem_dbg::MemSize::mem_size(
                        ef,
                        mem_dbg::SizeFlags::default()
                    ))
                    .sum::<usize>() as f32
                    / 1e9
            );

            eprintln!("Select pivots");
            let mut pivot_idxs = (0..1000)
                .map(|_| rand::random_range(0..num_vals))
                .collect_vec();
            pivot_idxs.sort_unstable();
            pivot_idxs.dedup();

            let mut pivots = vec![];
            {
                let mut i = 0;
                let mut pivot_idx_iter = pivot_idxs.iter();
                let mut next_pivot = *pivot_idx_iter.next().unwrap();
                'outer: for (_s, ef, _cn, _ce) in &child_results {
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
                pivots.push(0);
                pivots.push(u128::MAX);
                pivots.sort_unstable();
                pivots.dedup();
            }
            let mut efs_per_pivot = vec![vec![]; pivots.len()];
            let mut max = 0u128;
            eprintln!("Partitioning EFs");
            for (_s, ef, _cn, _ce) in child_results {
                if ef.len() == 0 {
                    continue;
                }
                let vals = ef.iter().collect_vec();
                max = max.max(*vals.last().unwrap());
                let splits = pivots
                    .iter()
                    .map(|&p| vals.partition_point(|&x| x < p))
                    .collect_vec();
                for i in 0..pivots.len() - 1 {
                    efs_per_pivot[i].push(BareEf::from(&vals[splits[i]..splits[i + 1]]));
                }
            }
            eprintln!(
                "Total EF size after partitioning: {:.3} GB",
                efs_per_pivot
                    .iter()
                    .flatten()
                    .map(|ef| mem_dbg::MemSize::mem_size(ef, mem_dbg::SizeFlags::default()))
                    .sum::<usize>() as f32
                    / 1e9
            );

            eprintln!("Build an EF for each part");
            let part_efs: Vec<(u128, BareEf)> = efs_per_pivot
                .into_par_iter()
                .enumerate()
                .map_with(vec![], |vals: &mut Vec<u128>, (i, efs)| {
                    eprintln!(
                        "{i}: Merging {} EFs of total len {} total size {:.3}",
                        efs.len(),
                        efs.iter().map(|ef| ef.len()).sum::<usize>(),
                        efs.iter()
                            .map(|ef| mem_dbg::MemSize::mem_size(ef, mem_dbg::SizeFlags::default()))
                            .sum::<usize>() as f32
                            / 1e9,
                    );
                    for ef in efs {
                        vals.extend(ef.iter());
                    }
                    eprintln!(
                        "{i}: vals size: {:.3} GB",
                        std::mem::size_of_val(vals.as_slice()) as f32 / 1e9
                    );
                    vals.sort_unstable();
                    vals.dedup();
                    let min = *vals.first().unwrap_or(&0);
                    for x in &mut *vals {
                        *x -= min;
                    }
                    let out = BareEf::from(&mut *vals);
                    vals.clear();
                    eprintln!(
                        "{i}: output EF size: {:.3} GB",
                        mem_dbg::MemSize::mem_size(&out, mem_dbg::SizeFlags::default()) as f32
                            / 1e9
                    );
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

            let n = part_efs.iter().map(|ef| ef.1.len()).sum();
            eprintln!("Merge part EFs. Dedupped to {n} links");
            let mut ef_builder = sux::dict::elias_fano::EliasFanoBuilder::<u128>::new(n, max);
            for (min, part_ef) in part_efs {
                eprintln!(
                    "Extending by {} elems size {}",
                    part_ef.len(),
                    mem_dbg::MemSize::mem_size(&part_ef, mem_dbg::SizeFlags::default()) as f32
                        / 1e9
                );
                ef_builder.extend(part_ef.iter().map(|x| x + min));
            }
            ef_builder.build_with_dict()
        };

        eprintln!(
            "final EF size: {:.3} GB",
            mem_dbg::MemSize::mem_size(&ef_links, mem_dbg::SizeFlags::default()) as f32 / 1e9
        );

        // eprintln!(
        //     "Links: {:.3} GB",
        //     std::mem::size_of_val(links.as_slice()) as f32 / 1e9
        // );

        let State { t, sa: _sa, .. } = state;

        // use voracious_radix_sort::RadixSort;
        // stpd_samples.voracious_mt_sort(12);
        // // FIXME: THIS IS TERRIBLY SLOW FOR 12 BYTE DATA.
        // links.voracious_mt_sort(12);

        stpd_samples.dedup();
        // links.dedup();
        // Free the excess capacity.
        // links.shrink_to_fit();
        // eprintln!(
        //     "Links: {:.3} GB (deduped)",
        //     std::mem::size_of_val(links.as_slice()) as f32 / 1e9
        // );

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
                let l = Link::from_key(*l);
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

        stpd_samples.sort_by(|&a, &b| cmp_colex(&t.as_ref()[..=a], &t.as_ref()[..=b]).1);
        let stpd_pi: Vec<u64> = stpd_samples.iter().map(|&x| pi[x] as u64).collect();
        // stpd_samples.iter().take(10).for_each(|&i| {
        //     eprintln!(
        //         "{i:>3}: {} ({})",
        //         crate::print(&t.as_ref()[i.saturating_sub(30)..=i]),
        //         pi[i]
        //     );
        // });

        JumpIndex {
            t,
            stpd_samples,
            stpd_rmq: rmq::BlockRmq::build(stpd_pi.as_slice()),
            stpd_pi,
            ef_links,
            cdawg_nodes,
            cdawg_edges,
        }
    }

    pub fn stats(&self) -> JumpIndexStats {
        JumpIndexStats {
            num_sampled: 1 + self
                .stpd_samples
                .iter()
                .tuple_windows()
                .filter(|(a, b)| a != b)
                .count(),
            num_sources: 1 + self
                .ef_links
                .iter()
                .tuple_windows()
                .filter(|&(a, b)| Link::from_key(a).source() != Link::from_key(b).source())
                .count(),
            num_source_chars: 1 + self
                .ef_links
                .iter()
                .tuple_windows()
                .filter(|&(a, b)| Link::from_key(a).source_c() != Link::from_key(b).source_c())
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
            "stpd samples {:.3} GB",
            std::mem::size_of_val(self.stpd_samples.as_slice()) as f32 / 1e9
        );
        eprintln!(
            "stpd pi      {:.3} GB",
            std::mem::size_of_val(self.stpd_pi.as_slice()) as f32 / 1e9
        );
        // eprintln!(
        //     "stpd rmq     {:.3} GB",
        //     std::mem::size_of_val(self.stpd_rmq.as_slice()) as f32 / 1e9
        // );
        eprintln!(
            "jump index   {:.3} GB",
            mem_dbg::MemSize::mem_size(&self.ef_links, mem_dbg::SizeFlags::default()) as f32 / 1e9
        );
    }

    /// Returns leftmost text position where the pattern matches.
    pub fn map_stpd(&self, pattern: &[u8]) -> Option<usize> {
        let mut pos = 0;
        for (i, &c) in pattern.iter().enumerate() {
            if self.t.as_ref()[pos] == c {
                pos += 1;
                continue;
            }

            // TODO: Use binary search function that reuses LCP.
            let idx1 = self
                .stpd_samples
                .binary_search_by(|&sample_pos| {
                    match cmp_colex(&self.t.as_ref()[..=sample_pos], &pattern[..=i]).1 {
                        std::cmp::Ordering::Equal => Greater,
                        x => x,
                    }
                })
                .unwrap_err();
            let idx2 = self
                .stpd_samples
                .binary_search_by(|&sample_pos| {
                    match cmp_colex(&self.t.as_ref()[..=sample_pos], &pattern[..=i]).1 {
                        std::cmp::Ordering::Equal => Less,
                        x => x,
                    }
                })
                .unwrap_err();
            if idx1 == idx2 {
                // eprintln!("idx: {idx1}..{idx2}");
                return None;
            }
            let (_val, idx) = self.stpd_rmq.query(self.stpd_pi.as_slice(), idx1, idx2 - 1);
            // eprintln!("idx: {idx1}..{idx2} => {idx} val={_val}");
            pos = self.stpd_samples[idx] + 1;

            // eprintln!("pos: {pos}");
        }
        // eprintln!("end pos {pos}");
        Some(pos - pattern.len())
    }

    /// Returns leftmost text position where the pattern matches.
    pub fn map_jump(&self, pattern: &[u8]) -> Option<usize> {
        let mut pos = 0;
        for (i, &c) in pattern.iter().enumerate() {
            if self.t.as_ref()[pos] == c {
                pos += 1;
                continue;
            }

            let (_idx, link) = self.ef_links.succ(&Link::new(pos, c, i, 0).key()).unwrap();
            let link = Link::from_key(link);
            // eprintln!("pos {pos} link {link:?}");
            if link.source() == pos && link.c() as u8 == c {
                pos = link.target() + 1;
            } else {
                return None;
            }
        }
        Some(pos - pattern.len())
    }
}

fn co_lcp(a: &[u8], b: &[u8]) -> usize {
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
    fn test() {
        let mut t1 = std::time::Duration::default();
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
            let t = relative(len, 4, repeats, r).1;
            eprintln!("text: {}", print(&t));

            eprintln!("building for {len}x{repeats} at {r}..");
            let ji = JumpIndex::new(&t);

            // find a bunch of random substrings
            eprintln!("querying..");
            for len in 0..=len.min(1000) {
                let pos = rand::random_range(0..=t.len() - len);
                let pattern = &t[pos..pos + len];

                eprintln!("pattern: {}", print(pattern));

                let s = std::time::Instant::now();
                let p1 = ji.map_stpd(pattern);
                t1 += s.elapsed();
                let s = std::time::Instant::now();
                let p2 = ji.map_jump(pattern);
                t2 += s.elapsed();
                eprintln!("p1: {p1:?}");
                eprintln!("p2: {p2:?}");
                let p1 = p1.unwrap();
                let p2 = p2.unwrap();

                assert_eq!(&t[p1..p1 + len], pattern);
                assert_eq!(&t[p2..p2 + len], pattern);
                assert_eq!(p1, p2);
            }
        }
        eprintln!("STPD: {t1:?}");
        eprintln!("JI:   {t2:?}");
    }
}
