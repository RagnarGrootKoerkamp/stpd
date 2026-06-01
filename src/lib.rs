#![feature(gen_blocks, bstr, vec_from_fn)]
use std::{cmp::{Ordering, Reverse}, collections::{BTreeSet, HashMap, HashSet, hash_map::Entry}};
use itertools::Itertools;
use rand::{rng, seq::SliceRandom};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

pub mod strings;
pub mod test;
pub mod stpd;

pub type T = Vec<u8>;
pub type SA = Vec<usize>;
pub type LCP = Vec<usize>;

pub fn print(t: &[u8]) -> String {
    if t[0] <= 4 {
        String::from_utf8(t.iter().map(|c| b'0' + c).collect_vec()).unwrap()
    } else {
        String::from_utf8(t.to_vec()).unwrap()
    }
}

pub fn sa(t: &T) -> SA {
    libsais::SuffixArrayConstruction::for_text(t.as_slice())
        .in_owned_buffer32()
        .multi_threaded(libsais::ThreadCount::openmp_default())
        .run()
        .unwrap()
        .suffix_array()
        .iter()
        .map(|&x| x as usize)
        .collect()
}

fn co_sa(t: &T) -> Vec<usize> {
    let tr = t.iter().rev().copied().collect_vec();
    let co_sa = sa(&tr);
    co_sa.into_iter().map(|x| t.len() - 1 - x).collect_vec()
}

pub fn sa_and_lcp(t: &T) -> (SA, LCP) {
    let sa_builder = libsais::SuffixArrayConstruction::for_text(t.as_slice())
        .in_owned_buffer32()
        .multi_threaded(libsais::ThreadCount::openmp_default())
        .run()
        .unwrap();
    let sa: Vec<usize> = sa_builder.suffix_array()
        .iter()
        .map(|&x| x as usize)
        .collect();

    let n = t.len();
    let lcp_raw = sa_builder.plcp_construction()
    .multi_threaded(libsais::ThreadCount::openmp_default())
    .run()
    .unwrap()
    .lcp_construction()
    .multi_threaded(libsais::ThreadCount::openmp_default())
    .run()
    .unwrap();
    let (_, lcp_raw, _, _) = lcp_raw.into_parts();
    // Drop the sentinel at index 0; the remaining n-1 values correspond to lcp[0..n-1].
    let mut result: LCP = lcp_raw[1..].iter().map(|&x| x as usize).collect();
    // Append the circular wrap: LCP(sa[n-1], sa[0]).
    let (a, b) = (sa[n - 1], sa[0]);
    let mut l = 0;
    while a + l < n && b + l < n && t[a + l] == t[b + l] {
        l += 1;
    }
    result.push(l);

    (sa, result)
}

pub fn bwt(t: &T, sa: &SA) -> T {
    sa.par_iter()
        .map(|&i| t[if i == 0 { t.len() - 1 } else { i - 1 }])
        .collect()
}

/// Number of BWT runs.
pub fn r(bwt: &T) -> usize {
    1 + bwt
        .iter()
        .circular_tuple_windows()
        .map(|(l, r)| if l != r { 1 } else { 0 })
        .sum::<usize>()
}

/// delta, and maximizing k.
pub fn delta(t: &T) -> (f32, usize) {
    let mut max = (0.0, 0);
    for k in 1..=t.len() {
        let mut kmers = HashSet::new();
        for kmer in t.windows(k) {
            kmers.insert(kmer);
        }
        let delta_k = (kmers.len() as f32 / k as f32, k);
        if delta_k > max {
            max = delta_k;
        }
    }
    max
}

/// Iterate right maximal strings α, corresponding to paths to explicit suffix tree nodes.
pub fn tree_nodes<'t>(t: &'t T, sa: &SA, lcp: &LCP) -> impl Iterator<Item = &'t [u8]> {
    gen {
        let mut depths = vec![0];
        for i in 0..t.len() {
            let l = lcp[i];
            while *depths.last().unwrap() > l {
                let l2 = depths.pop().unwrap();
                yield &t[sa[i]..sa[i]+l2];
            }
            if *depths.last().unwrap() < l {
                depths.push(l);
            }
        }
        yield b"";
    }
}

/// Iterate right maximal extensions αX, corresponding to edges going out of explicit suffix tree nodes.
pub fn tree_edges<'t>(t: &'t T, sa: &SA, lcp: &LCP) -> impl Iterator<Item = &'t [u8]> {
    gen {
        let mut depths = vec![];
        for i in 0..t.len() {
            let l = lcp[i];
            // eprintln!("{i} {l} {} {:?}", sa[i], &t[sa[i]..]);
            while let Some(l2) = depths.last() && *l2 > l {
                let l2 = depths.pop().unwrap();
                yield &t[sa[i]..sa[i]+l2+1];
            }
            // Otherwise the suffix is just an internal node anyway.
            if sa[i] + l + 1 <=  t.len(){
                yield &t[sa[i]..sa[i]+l+1];
            }
            depths.push(l);
        }
    }
}

/// Number of leaves in the Weiner-link-tree on explicit suffix tree nodes.
pub fn w(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut internal = HashSet::new();
    let mut num_nodes = 0;
    for n in tree_nodes(t,sa,lcp) {
        num_nodes += 1;
        if !n.is_empty() {
            internal.insert(&n[1..]);
        }
    }
    num_nodes - internal.len()
}

/// Size of smallest suffixient set.
pub fn chi(t: &T, sa: &SA, lcp: &LCP, print: bool) -> usize {
    let mut edges = tree_edges(t, sa, lcp).collect_vec();
    edges.sort_by_key(|e| Reverse(e.len()));

    let mut chi = 0;
    let mut covered = HashSet::new();
    covered.insert([].as_slice());
    for mut e in edges {
        if !covered.contains(e) {
            chi += 1;

            if print {
                eprintln!("{}", crate::print(e));
            }

            while !covered.contains(e) {
                covered.insert(e);
                e = &e[1..];
            }
        }
    }
    chi
}

// Returns a map alphaX -> alphaX..., where each right-maximal extension is mapped to its node.
pub fn tree<'t>(t: &'t T, sa: &SA, lcp: &LCP) -> HashMap<&'t [u8], &'t[u8]> {
    let mut nodes = HashSet::new();
    for n in tree_nodes(t, sa, lcp) {
        nodes.insert(n);
    }
    let mut tree = HashMap::new();
    for e in tree_edges(t, sa, lcp) {
        // Iterate up until a node is found.
        let mut p = e;
        while !nodes.contains(p.split_last().unwrap().1) {
            p.split_off_last().unwrap();
        }
        tree.insert(p, e);
    }
    tree
}

/// Size of smallest suffixient set, with path decomposition.
pub fn chi_pd(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut edges = tree_edges(t, sa, lcp).collect_vec();
    edges.sort_by_key(|e| Reverse(e.len()));

    let tree = tree(t, sa, lcp);

    let mut chi = 0;
    let mut covered = HashSet::new();
    let mut path_covered = HashSet::new();
    covered.insert([].as_slice());
    for mut e in edges {

        if !covered.contains(e) && !path_covered.contains(e) {
            chi += 1;

            while !covered.contains(e) {
                covered.insert(e);

                // Also cover right-extensions of e, after walking down the tree.
                {
                    // Position of start of e.
                    let pos = unsafe { e.as_ptr().offset_from(t.as_ptr()) as usize };
                    let mut e = &e[..e.len()-1];
                    while pos + e.len()+1 < t.len() {
                        // right-extend by 1 char to get to a right-maximal-extension.
                        e = &t[pos..pos+e.len()+1];
                        path_covered.insert(e);
                        if let Some(e2) = tree.get(e) {
                            e = e2;
                        } else {
                            break;
                        }

                    }
                }


                e = &e[1..];
            }
        }
    }
    chi
}

/// Size of smallest suffixient set, with path decomposition.
pub fn chi_pd2(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut edges = tree_edges(t, sa, lcp).collect_vec();
    edges.sort_by_key(|e| Reverse(e.len()));

    let tree = tree(t, sa, lcp);

    let mut chi = 0;
    let mut covered = HashSet::new();
    let mut path_covered = HashSet::new();
    covered.insert([].as_slice());
    for mut e in edges {

        if !covered.contains(e) && !path_covered.contains(e) {
            chi += 1;

            while e != b"" {
                covered.insert(e);

                // Also cover right-extensions of e, after walking down the tree.
                {
                    // Position of start of e.
                    let pos = unsafe { e.as_ptr().offset_from(t.as_ptr()) as usize };
                    let mut e = &e[..e.len()-1];
                    while pos + e.len()+1 < t.len() {
                        // right-extend by 1 char to get to a right-maximal-extension.
                        e = &t[pos..pos+e.len()+1];
                        path_covered.insert(e);
                        if let Some(e2) = tree.get(e) {
                            e = e2;
                        } else {
                            break;
                        }

                    }
                }


                e = &e[1..];
            }
        }
    }
    chi
}

pub fn stpd(t: &T, _sa: &SA, _lcp: &LCP, perm: &Vec<usize>) -> usize {
    let mut iperm = vec![0; t.len()];
    for i in 0..t.len(){
        iperm[perm[i]] = i;
    }
    // Map substrings to (sampled STPD pos, end)
    let mut seen = HashMap::<&[u8], (usize,usize)>::new();

    for i in 0..=t.len() {
        seen.insert(&t[0..i], (0,i));
    }

    // Map sampled STPD pos to (min length, max length, children, parents, suffix links, forward links)
    let mut sampled = HashMap::<usize, (usize, usize, HashSet<usize>, Vec<(usize,usize)>, Vec<(usize, usize)>)>::new();
    let mut fwd_links = HashMap::<usize, HashMap<u8, HashSet<usize>>>::new();

    // let mut branch_points = HashSet::new();
    // let mut parents = HashSet::new();
    // let mut suffix_links = HashSet::new();

    // eprintln!("T: {}", crate::print(t));
    for len in 1..=t.len(){
        for &pos in &iperm {
            if pos < len-1 {continue;}
            let start = pos-len+1;
            // eprintln!("t[{:?}] = {}", start..=pos, crate::print(&t[start..=pos]));
            let e = &t[start..=pos];
            if seen.contains_key(e) {
                continue;
            }
            // eprintln!("pos {pos:>2}: First time seeing {}", crate::print(e));

            // Add STPD element.
            let e = sampled.entry(pos).or_insert((len, len, Default::default(), vec![], Default::default()));
            // Update max len at current sample.
            e.1 = len;

            // let mut last_sl = usize::MAX;
            // let mut cnt = 0;
            for end in pos+1..=t.len(){
                // Mark seqs starting here as seen.
                seen.insert(&t[start..end], (pos, end));

                // Add suffix links.
                // if let Some(&(parent_pos, _parent_end)) = seen.get(&t[start+1..end]) {
                //     if parent_pos != last_sl && end-1 != pos {
                //         // Add suffix link to parent.
                //         e.4.push((end-1, parent_pos));
                //         cnt += 1;
                //         // Add suffix link to parent.
                //         // suffix_links.insert();
                //     }
                //         last_sl = parent_pos;
                // }
            }
            // assert!(cnt <= 2, "CNT: {cnt} at {pos}");

            if let Some((_parent_pos, parent_end)) = seen.get(&t[start..pos]) {
                // Add parent of this sample.
                // e.3.push((len-1, *parent_pos));
                // parents.insert((pos, *parent_pos));

                // Add this as child of previous sample.
                // sampled.get_mut(parent_pos).unwrap().2.insert(parent_end+1);
                // branch_points.insert(parent_end+1);

                // eprintln!("Parent {_parent_pos}..{parent_end}: {} => {pos}", t[pos]);
                // eprintln!("Add link from {parent_end} to {pos} for {}", t[pos] as char);
                fwd_links.entry(*parent_end).or_default().entry(t[pos]).or_default().insert(pos);
            }

        }
    }

    let mut sampled = sampled.into_iter().collect_vec();
    sampled.sort_unstable_by_key(|(pos, _)| *pos);
    // eprintln!("T: {}", crate::print(t));
    // for (pos, (min_len, max_len, children, parents, suffix_links)) in &sampled {
    //     eprintln!("{pos}: {min_len}..{max_len}  children: {:?}  parents: {:?}  suffix_links: {:?}", children, parents, suffix_links);
    // }
    // println!("{sampled:?}");
    eprint!("STPD samples:  {:>5}  | ", sampled.len());
    // eprint!("Branch points: {:>5}  | ", branch_points.len());
    // eprint!("Parents:       {:>5}  | ", parents.len());
    eprint!("Links 1:         {:>5}  | ", fwd_links.values().len());
    eprint!("Links 2:         {:>5}  | ", fwd_links.values().map(|m| m.len()).sum::<usize>());
    let l3 = fwd_links.values().map(|m|m.values().map(|x|x.len()).sum::<usize>()).sum:: <usize>();
    eprint!("Links 3:         {l3:>5}  | ");
    let l3max = fwd_links.values().map(|m|m.values().map(|x|x.len()).max().unwrap()).max().unwrap();
    eprint!("max:         {l3max:>5}  | ");
    eprintln!(" {:1.4}x", l3 as f32 / sampled.len() as f32);
    // eprintln!("Suffix links:  {:?}", suffix_links.len());
    sampled.len()
}

pub fn stpd_fast(t: &T, sa: &SA, bwt: &T, lcp: &LCP, pi: &Vec<usize>) -> usize {
    const PARALLEL_THRESHOLD: usize = 100_000;
    // eprintln!("T:   {}", crate::print(t));
    // eprintln!("sa:  {sa:?}");
    // eprintln!("lcp: {lcp:?}");
    // eprintln!("pi:  {pi:?}");
    // eprintln!("ppi: {permuted_pi:?}");
    struct State<'a> {
        t: &'a T,
        sa: &'a SA,
        lcp: &'a LCP,
        run_boundaries: BTreeSet<usize>,
        lcp_rmq: rmq_rust::BlockRmq<'a, 128>,
        #[allow(unused)]
        permuted_pi: &'a Vec<usize>,
        pi_rmq: rmq_rust::BlockRmq<'a, 128>,
    }

    impl<'a> State<'a> {
        /// Split an interval into its anchor, current LCP depth, and child intervals.
        fn split(&self, interval: std::ops::Range<usize>) -> Option<(usize, usize, Vec<std::ops::Range<usize>>)> {
            if interval.len() <= 1 {
                return None;
            }
            // if contained in a single BWT run, skip.
            // FIXME: We loose 1-2 samples with this. Find out exactly which.
            // NOTE: Understand why we can't move this into the loop below.
            if self.run_boundaries.range(interval.clone()).next().is_none() {
                return None;
            }
        
            let anchor_pos = self.pi_rmq.query(interval.start, interval.end - 1).1;
            let mut done_intervals = vec![];
            let mut wip_intervals = vec![interval.clone()];
            let lcp = self.lcp[self.lcp_rmq.query(interval.start, interval.end - 2).1];
            while let Some(interval) = wip_intervals.pop() {
                if interval.len() <= 1 {
                    done_intervals.push(interval);
                    continue;
                }
                let split_pos = self.lcp_rmq.query(interval.start, interval.end - 2).1 + 1;
                let new_lcp = self.lcp[split_pos - 1];
                if new_lcp > lcp {
                    done_intervals.push(interval);
                    continue;
                }
                assert!(new_lcp == lcp);
                wip_intervals.push(interval.start..split_pos);
                wip_intervals.push(split_pos..interval.end);
            }
            Some((anchor_pos, lcp, done_intervals))
        }

        /// Compute the sampled positions and forward links emitted at one suffix-tree node.
        fn node_output(&self, anchor_pos: usize, lcp: usize, done_intervals: &[std::ops::Range<usize>], sampled: &mut Vec<usize>, links: &mut Vec<Link>) {
            for x in done_intervals {
                if !x.contains(&anchor_pos) {
                    let secondary_anchor_pos = self.pi_rmq.query(x.start, x.end - 1).1;
                    let text_idx = self.sa[secondary_anchor_pos];
                    let target = text_idx + lcp;
                    if target < self.t.len() {
                        sampled.push(target);
                        let source = self.sa[anchor_pos] + lcp;
                        let c = self.t[target];
                        links.push(Link(source, c, target));
                    }
                }
            }
        }

        /// Sequential DFS used for each leaf task from the work queue.
        /// Results are accumulated directly into the caller-provided vecs.
        fn dfs(&self, interval: std::ops::Range<usize>, sampled: &mut Vec<usize>, links: &mut Vec<Link>) {
            let Some((anchor_pos, lcp, done_intervals)) = self.split(interval) else {
                return;
            };
            self.node_output(anchor_pos, lcp, &done_intervals, sampled, links);
            for x in done_intervals {
                self.dfs(x, sampled, links);
            }
        }

        /// Top-level sequential DFS: processes nodes normally until an interval falls below
        /// PARALLEL_THRESHOLD, then pushes it onto the work queue instead of recursing.
        fn collect_work(
            &self,
            interval: std::ops::Range<usize>,
            sampled: &mut Vec<usize>,
            links: &mut Vec<Link>,
            work_queue: &mut Vec<std::ops::Range<usize>>,
        ) {
            if interval.len() < PARALLEL_THRESHOLD {
                work_queue.push(interval);
                return;
            }
            let Some((anchor_pos, lcp, done_intervals)) = self.split(interval) else {
                return;
            };
            self.node_output(anchor_pos, lcp, &done_intervals, sampled, links);
            for x in done_intervals {
                self.collect_work(x, sampled, links, work_queue);
            }
        }
    }

    // indices where new runs start
    let run_boundaries = (0..t.len()-1).tuple_windows().filter(|(i, j)| bwt[*i] != bwt[*j]).map(|(i,_j)| i).collect_vec();
    let run_boundaries = std::collections::BTreeSet::from_iter(run_boundaries);

    let permuted_pi: Vec<usize> = sa.par_iter().map(|&i| pi[i]).collect();

    use rmq_rust::Rmq as _;
    let lcp_u64: Vec<u64> = lcp.iter().map(|&x| x as u64).collect();
    let ppi_u64: Vec<u64> = permuted_pi.iter().map(|&x| x as u64).collect();
    let state = State {
        t,
        sa,
        lcp,
        run_boundaries,
        lcp_rmq: rmq_rust::BlockRmq::build(&lcp_u64),
        pi_rmq: rmq_rust::BlockRmq::build(&ppi_u64),
        permuted_pi: &permuted_pi,
    };

    // Phase 1: sequential DFS collecting top-level results and a queue of leaf intervals.
    let mut sampled_vec = vec![];
    let mut links_vec = vec![];
    let mut work_queue = vec![];
    state.collect_work(0..t.len(), &mut sampled_vec, &mut links_vec, &mut work_queue);

    // Phase 2: process leaf intervals in parallel using Rayon's thread pool.
    // Each worker accumulates its results into a single pair of vecs.
    use rayon::prelude::*;
    let child_results: Vec<(Vec<usize>, Vec<Link>)> = work_queue
        .into_par_iter()
        .map(|interval| {
            let mut sampled = vec![];
            let mut links = vec![];
            state.dfs(interval, &mut sampled, &mut links);
            (sampled, links)
        })
        .collect();
    for (s, l) in child_results {
        sampled_vec.extend(s);
        links_vec.extend(l);
    }

    #[derive(Copy, PartialEq, PartialOrd, Clone)]
    struct Link(usize, u8, usize);
    impl voracious_radix_sort::Radixable<u128> for Link {
        type Key = u128;

        fn key(&self) -> Self::Key {
            ((self.0 as u128) << 72) | ((self.1 as u128) << 64) | (self.2 as u128)
        }
    }

    use voracious_radix_sort::RadixSort;
    sampled_vec.voracious_mt_sort(12);
    links_vec.voracious_mt_sort(12);

    let num_sampled = 1 + sampled_vec.iter().tuple_windows().filter(|(a, b)| a != b).count();
    let num_sources = 1 + links_vec.iter().tuple_windows().filter(|(a, b)| a.0 != b.0).count();
    let num_source_chars = 1 + links_vec.iter().tuple_windows().filter(|(a, b)| (a.0, a.1) != (b.0, b.1)).count();
    let num_links = 1 + links_vec.iter().tuple_windows().filter(|(a, b)| a != b).count();

    let c = 1000000.;
    let c = 1.;
    eprint!(" {:>5.2}  | ", num_sampled as f32 / c);
    eprint!(" {:>5.2}  | ", num_sources as f32 / c);
    eprint!(" {:>5.2}  | ", num_source_chars as f32 / c);
    eprintln!(" {:>5.2}  | ", num_links as f32 / c);
    num_sampled
}

pub fn stpd_pos_minus(t: &T, sa: &SA, bwt: &T, lcp: &LCP) -> usize {
    let perm = (0..t.len()).collect_vec();
    stpd_fast(t, sa, bwt, lcp, &perm)
}

pub fn stpd_pos_plus(t: &T, sa: &SA,  bwt: &T,lcp: &LCP) -> usize {
    let perm = (0..t.len()).rev().collect_vec();
    stpd_fast(t, sa, bwt, lcp, &perm)
}

pub fn stpd_lex_minus(t: &T, sa: &SA,  bwt: &T,lcp: &LCP) -> usize {
    let mut isa = vec![0; t.len()];
    for (i, &x) in sa.iter().enumerate(){
        isa[x] = i;
    }
    stpd_fast(t, sa, bwt, lcp, &isa)
}

pub fn stpd_lex_plus(t: &T, sa: &SA,  bwt: &T,lcp: &LCP) -> usize {
    let mut isa = vec![0; t.len()];
    for (i, &x) in sa.iter().enumerate(){
        isa[x] = t.len()-1-i;
    }
    stpd_fast(t, sa, bwt, lcp, &isa)
}

pub fn stpd_colex_minus(t: &T, sa: &SA,  bwt: &T,lcp: &LCP) -> usize {
    let co_sa = co_sa(t);
    let mut i_co_sa = vec![0; t.len()];
    for (i, &x) in co_sa.iter().enumerate(){
        i_co_sa[x] = i;
    }
    stpd_fast(t, sa, bwt, lcp, &i_co_sa)
}

pub fn stpd_colex_plus(t: &T, sa: &SA,  bwt: &T,lcp: &LCP) -> usize {
    let co_sa = co_sa(t);
    let mut i_co_sa = vec![0; t.len()];
    for (i, &x) in co_sa.iter().enumerate(){
        i_co_sa[x] = t.len()-1-i;
    }
    stpd_fast(t, sa, bwt, lcp, &i_co_sa)
}


pub fn stpd_rand(t: &T, sa: &SA, bwt: &T, lcp: &LCP) -> usize {
    let mut perm = (0..t.len()).collect_vec();
    perm.shuffle(&mut rng());
    stpd_fast(t, sa, bwt, lcp, &perm)
}

pub fn plcp(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut plcp = vec![0; t.len()];
    for i in 0..t.len(){
        plcp[sa[i]] = lcp[i];
    }
    let mut r = 1;
    for (&x, &y) in plcp.iter().tuple_windows() {
        if x != y+1 {
            r += 1;
        }
    }
    r
}
