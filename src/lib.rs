#![feature(gen_blocks, bstr, vec_from_fn, min_adt_const_params)]
use std::{cmp::{Reverse}, collections::{HashMap, HashSet}, hash::{Hash, Hasher}};
use itertools::Itertools;
use jump_index::{JumpIndexStats, Pi};
use lcp::{Lcp, PlainLcp};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};


use tikv_jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

pub mod strings;
pub mod test;
pub mod stpd;
pub mod jump_index;
mod rmq;
pub mod lcp;

pub type T = Vec<u8>;
pub type SaElem = u32;
pub type SA = Vec<u32>;
pub type LcpElem = u32;
pub type LCP = Vec<u32>;

pub fn print(t: &[u8]) -> String {
    if !t.is_empty() &&  t.iter().all(|x| *x <= 4) {
        String::from_utf8(t.iter().map(|c| b'0' + c).collect_vec()).unwrap()
    } else {
        String::from_utf8(t.to_vec()).unwrap()
    }
}

/// Shrink a Vec<i64> to a smaller type by iterating in reverse,
/// progressively freeing memory every million elements.
fn shrink_vec<T: Copy + Default>(mut source: Vec<i64>, convert: impl Fn(i64) -> T) -> Vec<T> {
    let len = source.len();
    let mut result = vec![T::default(); len];
    let shrink_interval = 1<<25;
    

    for i in (0..len).rev() {
        result[i] = convert(source[i]);
        if i % shrink_interval == 0 {
            source.truncate(i);
            source.shrink_to_fit();
        }
    }
    result
}

pub fn sa(t: &[u8]) -> SA {
    libsais::SuffixArrayConstruction::for_text(t)
        .in_owned_buffer64()
        .multi_threaded(libsais::ThreadCount::openmp_default())
        .run()
        .unwrap()
        .into_parts().0
        .into_iter()
        .map(|x| x as u32)
        .collect()
}

#[allow(unused)]
fn co_sa(t: &[u8]) -> SA {
    let tr = t.iter().rev().copied().collect_vec();
    let co_sa = sa(&tr);
    co_sa.into_iter().map(|x| t.len() as SaElem - 1 - x ).collect_vec()
}

type TheLcp = PlainLcp;

pub fn sa_and_lcp(t: &[u8]) -> (SA, TheLcp) {
    eprintln!("building sa..");
    let sa_builder = libsais::SuffixArrayConstruction::for_text(t)
        .in_owned_buffer64()
        .multi_threaded(libsais::ThreadCount::openmp_default())
        .run()
        .unwrap();
    eprintln!("SA:   {:.3} GB", std::mem::size_of_val(sa_builder.suffix_array()) as f32 / 1e9);

    eprintln!("building plcp..");
    let plcp_builder = sa_builder.plcp_construction()
    .multi_threaded(libsais::ThreadCount::openmp_default())
    .run()
    .unwrap();
    eprintln!("PLCP: {:.3} GB", std::mem::size_of_val(plcp_builder.plcp()) as f32 / 1e9);

    let (sa, plcp, _) = plcp_builder.into_parts();

    // Shrink types by progressively freeing allocations.
    eprintln!("shrinking..");
    let mut sa = shrink_vec(sa, |x| TryInto::<SaElem>::try_into(x).unwrap());
    eprintln!("SA:   {:.3} GB", std::mem::size_of_val(sa.as_slice()) as f32 / 1e9);
    let plcp = shrink_vec(plcp, |x|  TryInto::<LcpElem>::try_into(x).unwrap());
    eprintln!("PLCP: {:.3} GB", std::mem::size_of_val(plcp.as_slice()) as f32 / 1e9);

    // let plain_lcp = PlainLcp::new(&sa, &plcp);
    // eprintln!("PlainLCP: {:.3} GB", plain_lcp.space() as f32 / 1e9);

    // let compact_lcp = CompactLcp::new(&sa, &plcp);
    // eprintln!("CompactLCP: {:.3} GB", compact_lcp.space() as f32 / 1e9);

    // let ef_lcp = EfLcp::new(&sa, &plcp);
    // eprintln!("EfLCP: {:.3} GB", ef_lcp.space() as f32 / 1e9);

    let mut lcp = TheLcp::new(&sa, &plcp);
    eprintln!("{} LCP: {:.3} GB", std::any::type_name::<TheLcp>(), lcp.space() as f32 / 1e9);


    // Insert the empty suffix as smallest SA element.
    sa.insert(0, t.len() as SaElem);
    lcp.lcp.insert(0, 0);

    (sa, lcp)
}

pub fn sa_and_lcp_cached(t: &[u8]) -> (SA, TheLcp) {
    use std::collections::hash_map::DefaultHasher;
    use std::fs::{self, File};
    use std::io::{BufReader, BufWriter};
    use std::path::Path;

    // Hash the input text
    let mut hasher = DefaultHasher::new();
    t.hash(&mut hasher);
    let hash = hasher.finish();
    
    // Create cache directory if it doesn't exist
    let cache_dir = Path::new("_cache");
    if !cache_dir.exists() {
        fs::create_dir_all(cache_dir).expect("Failed to create _cache directory");
    }
    
    let cache_file = cache_dir.join(format!("{:x}.bin", hash));
    
    // Try to load from cache
    if cache_file.exists() {
        eprintln!("Loading from cache: {:?}", cache_file);
        let file = File::open(&cache_file).unwrap();
        let reader = BufReader::new(file);
        let (sa, lcp) = bincode::deserialize_from::<_, (SA, TheLcp)>(reader).unwrap();
        eprintln!("SA:   {:.3} GB", std::mem::size_of_val(sa.as_slice()) as f32 / 1e9);
        eprintln!("{}: {:.3} GB", std::any::type_name::<TheLcp>(), lcp.space() as f32 / 1e9);
        return (sa, lcp);
    }
    
    // Compute SA and LCP
    let result = sa_and_lcp(t);
    
    // Write to cache
    if t.len() > 100_000_000 {
        eprintln!("Writing to cache: {:?}", cache_file);
        let file = File::create(&cache_file).unwrap();
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &result).unwrap();
    }

    result
}

pub fn bwt(t: &[u8], sa: &SA) -> T {
    sa.par_iter()
        .map(|&i| t[if i == 0 { t.len() - 1 } else { i as usize - 1 }])
        .collect()
}

/// Number of BWT runs.
pub fn r(bwt: &[u8]) -> usize {
    1 + bwt
        .iter()
        .tuple_windows()
        .map(|(l, r)| if l != r { 1 } else { 0 })
        .sum::<usize>()
}

/// delta, and maximizing k.
pub fn delta(t: &[u8]) -> (f32, usize) {
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
pub fn tree_nodes<'t>(t: &'t T, sa: &SA, lcp: &impl Lcp) -> impl Iterator<Item = &'t [u8]> {
    gen {
        let mut depths = vec![0];
        for i in 0..t.len() {
            let l = lcp.get(sa, i);
            while *depths.last().unwrap() > l {
                let l2 = depths.pop().unwrap();
                yield &t[sa[i] as usize..sa[i] as usize+l2 as usize];
            }
            if *depths.last().unwrap() < l {
                depths.push(l);
            }
        }
        yield b"";
    }
}

/// Iterate right maximal extensions αX, corresponding to edges going out of explicit suffix tree nodes.
pub fn tree_edges<'t>(t: &'t T, sa: &SA, lcp: &impl Lcp) -> impl Iterator<Item = &'t [u8]> {
    gen {
        let mut depths = vec![];
        for i in 0..t.len() {
            let l = lcp.get(sa, i);
            // eprintln!("{i} {l} {} {:?}", sa[i], &t[sa[i]..]);
            while let Some(l2) = depths.last() && *l2 > l {
                let l2 = depths.pop().unwrap();
                let idx = sa[i] as usize;
                yield &t[idx..idx+l2 as usize+1];
            }
            // Otherwise the suffix is just an internal node anyway.
            if sa[i] as usize + l as usize + 1 <=  t.len(){
                let idx = sa[i] as usize;
                yield &t[idx..idx+l as usize+1];
            }
            depths.push(l);
        }
    }
}

/// Number of leaves in the Weiner-link-tree on explicit suffix tree nodes.
pub fn w(t: &T, sa: &SA, lcp: &impl Lcp) -> usize {
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
pub fn chi(t: &T, sa: &SA, lcp: &impl Lcp, print: bool) -> usize {
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
pub fn tree<'t>(t: &'t T, sa: &SA, lcp: &impl Lcp) -> HashMap<&'t [u8], &'t[u8]> {
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
pub fn chi_pd(t: &T, sa: &SA, lcp: &impl Lcp) -> usize {
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
pub fn chi_pd2(t: &T, sa: &SA, lcp: &impl Lcp) -> usize {
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

pub fn stpd(t: &T, _sa: &SA, _lcp: &impl Lcp, perm: &Vec<usize>) -> usize {
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

pub fn jump_index(t: &T) {
    let ji = jump_index::JumpIndex::<{Pi::LeftMost}>::new_dynamic(t);
    // let ji = jump_index::JumpIndex::<{Pi::LeftMost}>::new(t);

    // let JumpIndexStats { num_sampled, num_sources, num_source_chars, num_links, cdawg_nodes, cdawg_edges } = jump_index.stats();

    // #[allow(unused)]
    // let c = 1000000.;
    // let c = 1.;
    // eprint!(" {:>5.2}  | ", num_sampled as f32 / c);
    // eprint!(" {:>5.2}  | ", num_sources as f32 / c);
    // eprint!(" {:>5.2}  | ", num_source_chars as f32 / c);
    // eprint!(" {:>5.2}  | ", num_links as f32 / c);
    // eprint!(" {:>5.2}  | ", cdawg_nodes as f32 / c);
    // eprintln!(" {:>5.2}  | ", cdawg_edges as f32 / c);

    eprintln!("Testing mapping");
    ji.test_map();
    eprintln!("Testing MS");
    ji.test_ms();
    eprintln!("Bench MS");
    ji.bench_ms();
    eprintln!("Bench RLZ");
    ji.bench_rlz();

}


pub fn stpd_via_jump_index<const PI: Pi>(t: &T) -> usize {
    let jump_index = jump_index::JumpIndex::<PI>::new(t);

    let JumpIndexStats { num_sampled, num_sources, num_source_chars, num_links, cdawg_nodes, cdawg_edges } = jump_index.stats();

    #[allow(unused)]
    let c = 1000000.;
    let c = 1.;
    eprint!(" {:>5.2}  | ", num_sampled as f32 / c);
    eprint!(" {:>5.2}  | ", num_sources as f32 / c);
    eprint!(" {:>5.2}  | ", num_source_chars as f32 / c);
    eprint!(" {:>5.2}  | ", num_links as f32 / c);
    eprint!(" {:>5.2}  | ", cdawg_nodes as f32 / c);
    eprintln!(" {:>5.2}  | ", cdawg_edges as f32 / c);
    // jump_index.space();
    // jump_index.inspect_links();
    num_sampled
}

pub fn stpd_pos_minus(t: &T) -> usize {
    stpd_via_jump_index::<{Pi::LeftMost}>(t)
}

pub fn stpd_pos_plus(t: &T) -> usize {
    stpd_via_jump_index::<{Pi::RightMost}>(t)
}

pub fn stpd_lex_minus(_t: &T) -> usize {
    todo!();
}

pub fn stpd_lex_plus(_t: &T) -> usize {
    todo!();
}

pub fn stpd_colex_minus(_t: &T) -> usize {
    todo!();
}

pub fn stpd_colex_plus(_t: &T) -> usize {
    todo!();
}

pub fn plcp(t: &T, sa: &SA, lcp: &impl Lcp) -> usize {
    let mut plcp = vec![0; t.len()];
    for i in 0..t.len(){
        plcp[sa[i] as usize] = lcp.get(sa, i);
    }
    let mut r = 1;
    for (&x, &y) in plcp.iter().tuple_windows() {
        if x != y+1 {
            r += 1;
        }
    }
    r
}

pub fn longest_common_prefix(a: &[u8], b: &[u8]) -> usize {
    let l = a.len().min(b.len());
    for i in 0..l {
        if a[i] != b[i] {
            return i;
        }
    }
    l
}

pub fn gbs(ef: &impl mem_dbg::MemSize) -> f32 {
    mem_dbg::MemSize::mem_size(ef, mem_dbg::SizeFlags::default()) as f32 / 1e9
}
