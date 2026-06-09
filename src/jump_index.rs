use itertools::Itertools;
use mem_dbg::SizeFlags;
use std::{
    cmp::Ordering::{Greater, Less},
    collections::BTreeSet,
    marker::Sync,
};

use crate::{
    bwt,
    rmq::{self, Rmq},
    sa_and_lcp,
    stpd::cmp_colex,
    SaElem, LCP, SA, T,
};

#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Debug)]
pub struct Link {
    pub source: usize,
    pub c: char,
    pub lcp: usize,
    pub target: usize,
}

impl Link {
    fn key(&self) -> u128 {
        ((self.source as u128) << 64) | ((self.c as u128) << 48) | (self.lcp as u128)
    }
}

impl voracious_radix_sort::Radixable<u128> for Link {
    type Key = u128;
    fn key(&self) -> Self::Key {
        self.key()
    }
}

pub struct JumpIndex<TR: AsRef<T>> {
    pub t: TR,
    pub stpd_samples: Vec<usize>,
    pub stpd_pi: Vec<u64>,
    pub stpd_rmq: rmq::BlockRmq<u64, 64>,
    // TODO: Predecessor structure
    pub links: Vec<Link>,
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
        Self::new2(t, sa, bwt, lcp, &vec![])
    }
}
impl<TR: AsRef<T> + Sync> JumpIndex<TR> {
    pub fn new2<SAR: AsRef<SA> + Sync, LCPR: AsRef<LCP> + Sync>(
        t: TR,
        sa: SAR,
        bwt: &T,
        lcp: LCPR,
        pi: &SA,
    ) -> Self {
        const PARALLEL_THRESHOLD: usize = 100_000;

        struct State<'a, T, SA> {
            t: T,
            sa: &'a SA,
            lcp: &'a LCP,
            run_boundaries: Vec<SaElem>,
            lcp_rmq: rmq::BlockRmq<crate::LcpElem, 128>,
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
                let idx = self
                    .run_boundaries
                    .binary_search(&(interval.start as SaElem))
                    .unwrap_or_else(|x| x);
                let single_run = idx == self.run_boundaries.len()
                    || self.run_boundaries[idx] >= interval.end as SaElem - 1;
                if single_run {
                    return None;
                }

                let anchor_idx = self
                    .pi_rmq
                    .query(self.permuted_pi.as_ref(), interval.start, interval.end - 1)
                    .1;
                // let _anchor_pos = self.sa.as_ref()[anchor_idx];
                // eprintln!("anchor pos: {anchor_pos}");
                let mut done_intervals = vec![];
                let mut wip_intervals = vec![interval.clone()];
                let lcp = self.lcp[self
                    .lcp_rmq
                    .query(&self.lcp, interval.start, interval.end - 2)
                    .1];
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
                        .query(&self.lcp, interval.start, interval.end - 2)
                        .1
                        + 1;
                    let new_lcp = self.lcp[split_pos - 1];
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
                            .query(self.permuted_pi.as_ref(), x.start, x.end - 1)
                            .1;
                        let text_idx = self.sa.as_ref()[secondary_anchor_pos] as usize;
                        let target = text_idx + lcp as usize;
                        if target < self.t.as_ref().len() {
                            // NOTE: STPD samples are skipped for now.
                            // sampled.push(target);
                            let source = self.sa.as_ref()[anchor_pos] as usize + lcp as usize;
                            let c = self.t.as_ref()[target];
                            links.push(Link {
                                source,
                                c: c as char,
                                lcp: co_lcp(&self.t.as_ref()[..source], &self.t.as_ref()[..target]),
                                target,
                            });
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
                if interval.len() < PARALLEL_THRESHOLD {
                    work_queue.push(interval);
                    return;
                }
                let Some((single_run, anchor_pos, lcp, done_intervals)) = self.split(interval)
                else {
                    return;
                };
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

        let run_boundaries = (0..t.as_ref().len())
            .tuple_windows()
            .filter(|(i, j)| bwt[*i] != bwt[*j])
            .map(|(i, _j)| i as SaElem)
            .collect_vec();
        eprintln!(
            "Run boundaries: {:.3} GB",
            std::mem::size_of_val(run_boundaries.as_slice()) as f32 / 1e9
        );
        // let run_boundaries = BTreeSet::from_iter(run_boundaries);
        // let s = mem_dbg::MemSize::mem_size(&run_boundaries, SizeFlags::default());
        // eprintln!("Run boundaries set: {:.3} GB", s as f32 / 1e9);

        // empty indicates pi=identity and permuted_pi = sa.
        let permuted_pi: &SA = if pi.is_empty() {
            &sa.as_ref()
        } else {
            &sa.as_ref().par_iter().map(|&i| pi[i as usize]).collect()
        };

        // eprintln!();
        // eprintln!("sa:  {:?}", sa.as_ref());
        // eprintln!("lcp: {:?}", lcp.as_ref());
        // eprintln!("ppi: {:?}", permuted_pi);
        // eprintln!("run boundaries: {:?}", run_boundaries);

        use rmq::Rmq as _;
        let state = State {
            t,
            sa: sa.as_ref(),
            run_boundaries,
            lcp_rmq: rmq::BlockRmq::build(lcp.as_ref()),
            lcp: lcp.as_ref(),
            pi_rmq: rmq::BlockRmq::build(&permuted_pi),
            permuted_pi: &permuted_pi,
        };

        let mut stpd_samples = vec![];
        let mut links = vec![];
        let mut work_queue = vec![];
        let mut cdawg_nodes = 1; // the final node
        let mut cdawg_edges = 0;
        state.collect_work(
            0..state.t.as_ref().len(),
            &mut stpd_samples,
            &mut links,
            &mut work_queue,
            &mut cdawg_nodes,
            &mut cdawg_edges,
        );

        use rayon::prelude::*;
        let child_results: Vec<(Vec<usize>, Vec<Link>, usize, usize)> = work_queue
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
                (sampled, links, cdawg_nodes, cdawg_edges)
            })
            .collect();
        for (_s, l, cn, ce) in child_results {
            // stpd_samples.extend(s);
            links.extend(l);
            cdawg_nodes += cn;
            cdawg_edges += ce;
        }

        let State { t, sa: _sa, .. } = state;

        use voracious_radix_sort::RadixSort;
        stpd_samples.voracious_mt_sort(12);
        links.voracious_mt_sort(12);

        eprintln!(
            "Links: {:.3} GB",
            std::mem::size_of_val(links.as_slice()) as f32 / 1e9
        );
        stpd_samples.dedup();
        links.dedup();
        eprintln!(
            "Links: {:.3} GB (deduped)",
            std::mem::size_of_val(links.as_slice()) as f32 / 1e9
        );

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
            stpd_rmq: rmq::BlockRmq::build(&stpd_pi),
            stpd_pi,
            links,
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
                .links
                .iter()
                .tuple_windows()
                .filter(|(a, b)| a.source != b.source)
                .count(),
            num_source_chars: 1 + self
                .links
                .iter()
                .tuple_windows()
                .filter(|(a, b)| (a.source, a.c) != (b.source, b.c))
                .count(),
            num_links: 1 + self
                .links
                .iter()
                .tuple_windows()
                .filter(|(a, b)| a != b)
                .count(),
            cdawg_nodes: self.cdawg_nodes,
            cdawg_edges: self.cdawg_edges,
        }
    }

    pub fn inspect_links(&self) {
        eprintln!("Links: {}", self.links.len());
        if self.links.len() < 10000 {
            for i in 0..self.links.len() {
                eprintln!("{i:>8} {:?}", self.links[i]);
            }
        } else {
            for i in (0..self.links.len()).step_by(100000) {
                for i in i..i + 100 {
                    eprintln!("{i:>8} {:?}", self.links[i]);
                }
                eprintln!("---");
            }
        }
        eprintln!("---");

        eprintln!("Samples: {}", self.stpd_samples.len());
        if self.stpd_samples.len() < 30000 {
            for i in 0..self.stpd_samples.len() {
                eprintln!("{i:>8} {:>7}", self.stpd_samples[i]);
            }
        } else {
            for i in (0..self.stpd_samples.len()).step_by(100000) {
                for i in i..i + 100 {
                    eprintln!("{i:>8} {:>7}", self.stpd_samples[i]);
                }
                eprintln!("---");
            }
        }
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
            std::mem::size_of_val(self.links.as_slice()) as f32 / 1e9
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
            let (_val, idx) = self.stpd_rmq.query(&self.stpd_pi, idx1, idx2 - 1);
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

            let link_idx = self
                .links
                .binary_search_by(|link| {
                    link.key().cmp(
                        &Link {
                            source: pos,
                            c: c as char,
                            lcp: i,
                            target: 0,
                        }
                        .key(),
                    )
                })
                .map_or_else(|e| e, |v| v);
            let link = self.links[link_idx];
            // eprintln!("pos {pos} link {link:?}");
            if link.source == pos && link.c as u8 == c {
                pos = link.target + 1;
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
