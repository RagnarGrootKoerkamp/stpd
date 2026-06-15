#![allow(unused)]
#![feature(bstr)]

use std::bstr::ByteStr;

use text_indexing::strings::*;
use text_indexing::{test::Test, *};

fn header() {
    return;
    let name = "name";
    let n = "n";
    let r = "r";
    let delta = "δ";
    let delta_k = "δₖ";
    let delta_lg = "δlg(n/δ)";
    let w = "W";
    let chi = "χ";
    let chi_pd = "χ pd";
    let chi_pd2 = "χ set";
    let nodes = "nodes";
    let edges = "edges";
    let avg_node_depth = "n.dp";
    let avg_edge_depth = "e.dp";
    let inv_avg_node_depth = "in.dp";
    let inv_avg_edge_depth = "ie.dp";
    let normalized_tree_size = "ntsz";
    let stpd_pos_minus = "pos-";
    let stpd_pos_plus = "pos+";
    let stpd_lex_minus = "lex-";
    let stpd_lex_plus = "lex+";
    let stpd_colex_minus = "clex-";
    let stpd_colex_plus = "clex+";
    let stpd_rand = "rand";
    let plcp = "plcp";
    // {nodes:>5} {edges:>5}  \
    // {avg_node_depth:>4} {avg_edge_depth:>4} {inv_avg_node_depth:>5} {inv_avg_edge_depth:>5}  {normalized_tree_size:>4}  \
    eprintln!("{name:>40}  {n:>4}  \
{r:>4} {delta:>3} {delta_k:>3} {w:>5} {chi:>5} {chi_pd:>5} {chi_pd2:>5}  {delta_lg:>8} \
{stpd_pos_minus:>5} {stpd_pos_plus:>5} {stpd_lex_minus:>5} {stpd_lex_plus:>5} {stpd_colex_minus:>5} {stpd_colex_plus:>5} {stpd_rand:>5}");
}

fn stats((name, t): &(String, T), print: bool) {
    if print {
        eprintln!("T: {}", crate::print(t));
    }
    let n = t.len();
    let (sa, lcp) = &sa_and_lcp_cached(t);

    // let mut nodes = 0;
    // let mut node_depth = 0;
    // for n in tree_nodes(t, sa, lcp) {
    //     nodes += 1;
    //     node_depth += n.len();
    // }
    // let avg_node_depth = node_depth / nodes;
    // let inv_avg_node_depth = n * nodes / node_depth;

    // let mut edges = 0;
    // let mut edge_depth = 0;
    // for e in tree_edges(t, sa, lcp) {
    //     edges += 1;
    //     edge_depth += e.len();
    // }
    // let avg_edge_depth = edge_depth / edges;
    // let inv_avg_edge_depth = n * edges / edge_depth;

    // let mut tree_size = n * (n - 1) / 2;
    // for n in tree_nodes(t, sa, lcp) {
    //     tree_size += n.len() + 1;
    // }
    // for e in tree_edges(t, sa, lcp) {
    //     tree_size -= e.len();
    // }
    // let normalized_tree_size = tree_size / n;

    let bwt = &bwt(t, sa);
    // eprintln!("t:   {}", crate::print(t));
    // eprintln!("bwt: {}", crate::print(bwt));
    eprintln!(
        "t:   {:.3} GB",
        std::mem::size_of_val(t.as_slice()) as f32 / 1_000_000_000.
    );
    eprintln!(
        "bwt: {:.3} GB",
        std::mem::size_of_val(bwt.as_slice()) as f32 / 1_000_000_000.
    );
    let r = r(bwt);

    // // slow
    // let (delta, delta_k) = delta(t);
    // let delta_lg = (delta * (n as f32 / delta).log2()) as usize;
    // let delta = delta as usize;
    // let w = w(t, sa, lcp);
    // let chi = chi(t, sa, lcp, print && false);
    // let chi_pd = chi_pd(t, sa, lcp);
    // let chi_pd2 = chi_pd2(t, sa, lcp);
    // // slow

    let c = 1.;
    let n = n as f32 / c;
    let r = r as f32 / c;
    eprint!("| {name:<30} | {n:>6.2} | {r:>6.2} | pos-   | ");
    eprintln!();
    let stpd_pos_minus = stpd_pos_minus(t, sa, bwt, lcp);
    return;
    // eprint!("| {name} | {n:>6.2} | {r:>6.2} | pos+   | ");
    // let stpd_pos_plus = stpd_pos_plus(t, sa, bwt, lcp);
    eprint!("| {:<30} | {:>6} | {:>6} | lex-   | ", "", "", "");
    eprintln!();
    let stpd_lex_minus = stpd_lex_minus(t, sa, bwt, lcp);
    eprint!("| {:<30} | {:>6} | {:>6} | colex- | ", "", "", "");
    eprintln!();
    let stpd_colex_minus = stpd_colex_minus(t, sa, bwt, lcp);

    // // let plcp = plcp(t, sa, lcp);

    // // {nodes:>5} {edges:>5}  \
    // // {avg_node_depth:>4} {avg_edge_depth:>4} {inv_avg_node_depth:>5} {inv_avg_edge_depth:>5}  {normalized_tree_size:>4}  \
    // // eprintln!("{name:>40}  {n:>4}  \
    // // {r:>4} {delta:>3} {delta_k:>3} {w:>5} {chi:>5} {chi_pd:>5} {chi_pd2:>5}  {delta_lg:>8} \
    // // {stpd_pos_minus:>5} {stpd_pos_plus:>5} {stpd_lex_minus:>5} {stpd_lex_plus:>5} {stpd_colex_minus:>5} {stpd_colex_plus:>5} {stpd_rand:>5}");
}

fn stpd() {
    let len = 1_000_000;
    let r = 0.001;
    let muts = len as f32 * r;
    let start = std::time::Instant::now();
    let t = &relative(len, 4, 100, r).1;
    let duration = start.elapsed();
    log::error!("Gen took {duration:?}");

    let mut stpd = stpd::Stpd::new(b"");
    let mut a = 1;
    for (i, t) in t.chunks(len).enumerate() {
        let start = std::time::Instant::now();
        stpd.push(t);
        let duration = start.elapsed();
        let mbps = (t.len() as f64) / (duration.as_secs_f64() * 1_000_000.0);
        let added = stpd.num_anchors() - a;
        a = stpd.num_anchors();

        log::error!(
            "{i:>3}: {:.2} seconds ({:.2} Mbp/s)  {added:>6} new anchors {:4.1}/mut",
            duration.as_secs_f64(),
            mbps,
            added as f32 / muts,
        );
    }

    // RopeBWT: 65h for 320 copies of 3.2Gbp =>  4.2 Mbp/s  many threads (?)
    // us:      16s for 100 copies of 1  Mbp =>  6.2 Mbp/s  1 thread
    // us:     410s for 100 copies of 10 Mbp =>  2.4 Mbp/s  1 thread
    // now:
    // us:      10s for 100 copies of 1  Mbp => 10   Mbp/s  1 thread
    // us:     330s for 100 copies of 10 Mbp =>  3.0 Mbp/s  1 thread
}

fn stpd_human() {
    let mut stpd = stpd::Stpd::new(b"");
    let mut reader = needletail::parse_fastx_file("human-genome.fa").unwrap();
    while let Some(record) = reader.next() {
        // Convert to ascii ABCD
        let mut seq = record.unwrap().seq().to_vec();
        for b in seq.iter_mut() {
            *b = b'A' + ((*b >> 1) & 3);
        }

        let start = std::time::Instant::now();
        stpd.push(&seq);
        let duration = start.elapsed();
        let mbps = (seq.len() as f64) / (duration.as_secs_f64() * 1_000_000.0);
        log::error!(
            "{} MB {:.2} seconds ({:.2} Mbp/s)",
            seq.len() / 1_000_000,
            duration.as_secs_f64(),
            mbps
        );
    }

    // RopeBWT: 65h for 320 copies of 3.2Gbp =>  4.2 Mbp/s  many threads (?)
    // us:      16s for 100 copies of 1  Mbp =>  6.2 Mbp/s  1 thread
    // us:     410s for 100 copies of 10 Mbp =>  2.4 Mbp/s  1 thread
    // now:
    // us:      10s for 100 copies of 1  Mbp => 10   Mbp/s  1 thread
    // us:     330s for 100 copies of 10 Mbp =>  3.0 Mbp/s  1 thread
}

fn to_dna(mut t: Vec<u8>) -> Vec<u8> {
    for x in &mut t {
        *x = (*x >> 1) & 3;
    }
    t
}

fn pizzachili(filter: Option<&str>) -> Vec<(String, T)> {
    let dir = "/home/philae/git/eth/data/pizzachili/repetitive";
    let mut entries: Vec<_> = std::fs::read_dir(dir)
        .expect("failed to read pizzachili dir")
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.file_type()
                .map(|t| t.is_file() || t.is_symlink())
                .unwrap_or(false)
        })
        .filter(|e| {
            filter
                .map(|f| e.file_name().to_string_lossy() == f)
                .unwrap_or(true)
        })
        .collect();
    entries.sort_by_key(|e| e.file_name());
    entries
        .into_iter()
        .map(|e| {
            let name = e.file_name().to_string_lossy().into_owned();
            let data = std::fs::read(e.path()).expect("failed to read file");
            let data = to_dna(data);
            (name, data)
        })
        .collect()
}

fn ragc() -> Vec<(String, T)> {
    let path = "/home/philae/git/eth/data/hprcv2.agc";
    use ragc_core::{Decompressor, DecompressorConfig};

    // Open an archive
    let config = DecompressorConfig::default();
    let mut decompressor = Decompressor::open(path, config).unwrap();

    // List available samples
    let samples = decompressor.list_samples();
    println!("Found {} samples", samples.len());
    eprintln!("samples: {samples:?}");

    // Extract a sample
    let sample = &samples[1];
    let contigs = decompressor.list_contigs(sample).unwrap();
    eprintln!("contigs {contigs:?}");
    for contig in &contigs {
        let contig = decompressor.get_contig(sample, contig).unwrap();
    }
    // for (name, sequence) in contigs {
    //     println!(">{}", name);
    //     // sequence is Vec<u8> with numeric encoding (A=0, C=1, G=2, T=3)
    // }
    panic!();
}

fn main() {
    env_logger::Builder::from_default_env()
        .format_timestamp_micros()
        .init();

    // return stpd_human();
    // return stpd();

    // newtest();

    // return;

    // stats(&terminate(fib(6)), true);
    // return;

    let dataset = std::env::args().nth(1);

    header();
    let repeated = relative(200, 4, 20, 0.05);
    let texts = [
        // ragc(),
        pizzachili(dataset.as_deref()),
        // vec![("manual".to_string(), b"AGAGCGAGAGCGCGC#".to_vec())],
        // variants(fib(15)),
        vec![
            // relative(100, 4, 200, 0.01),
            // random(3_200_000_00, 4),
            // thue_morse(10),
            // random(100, 4),
            // relative(500, 4, 2, 0.01),
            // relative(250, 4, 4, 0.01),
            // relative(100, 4, 10, 0.01),
            // relative(50, 4, 20, 0.01),
            // relative(25, 4, 40, 0.01),
            // relative(250, 4, 5, 0.01),
            // relative(250, 4, 10, 0.01),
            // relative(250, 4, 20, 0.01),
            // relative(100, 4, 5, 0.01),
            // relative(100, 4, 10, 0.01),
            // relative(100, 4, 20, 0.01),
            // relative(100, 4, 40, 0.01),
            // relative(10, 4, 1, 0.01),
            // relative(10, 4, 2, 0.01),
            // relative(10, 4, 4, 0.01),
            // relative(10, 4, 8, 0.01),
            // relative(10000, 2, 1, 0.001),
            // relative(10000, 2, 2, 0.001),
            // relative(10000, 2, 128, 0.001),
            // relative(10000, 2, 1024, 0.001),
            // relative(10000, 4, 1, 0.001),
            // relative(10000, 4, 2, 0.001),
            // relative(10000, 4, 128, 0.001),
            // relative(10000, 4, 1024, 0.001),
            // relative(100000, 4, 1, 0.00),
            // relative(100000, 4, 2, 0.00),
            // relative(100000, 4, 3, 0.00),
            // relative(100000, 4, 4, 0.00),
            // relative(100000, 4, 1, 0.001),
            // relative(100000, 4, 2, 0.001),
            // relative(100000, 4, 3, 0.001),
            // relative(100000, 4, 4, 0.001),
            // relative(100000, 4, 1, 0.01),
            // relative(100000, 4, 2, 0.01),
            // relative(100000, 4, 3, 0.01),
            // relative(100000, 4, 4, 0.01),
            // relative(1000000, 4, 1, 0.001),
            // relative(1000000, 4, 2, 0.001),
            // relative(1000000, 4, 3, 0.001),
            // relative(1000000, 4, 4, 0.001),
            // relative(1000000, 4, 1, 0.01),
            // relative(1000000, 4, 2, 0.01),
            // relative(1000000, 4, 3, 0.01),
            // relative(1000000, 4, 4, 0.01),
            // relative(50, 4, 8, 0.02),
            // relative(50, 4, 20, 0.05),
            // relative(25, 4, 40, 0.05),
            // nice
            // relative(50, 4, 8, 0.02),
            // relative(25, 2, 100, 0.02),
        ],
        // vec![
        //     repeated.clone(),
        //     u8_minimizers(repeated.clone(), 1),
        //     u8_minimizers(repeated.clone(), 2),
        //     u8_minimizers(repeated.clone(), 4),
        //     u8_minimizers(repeated.clone(), 8),
        //     u8_minimizers(repeated.clone(), 16),
        //     u8_minimizers(repeated.clone(), 32),
        // ],
    ]
    .concat();
    for t in texts {
        // stats(&t, false);
        jump_index(&t.1);
    }

    // stats(&random(1000, 2), true);
}

fn newtest() {
    // let t = b"0100101001001010010100100101001001$".as_slice();
    let t = b"01001010$".as_slice();
    // let t = random(100, 2).1;
    Test::new(&t);
    eprintln!("{}", print(&t));
}
