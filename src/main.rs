#![allow(unused)]

use text_indexing::strings::*;
use text_indexing::{test::Test, *};

fn header() {
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
    let sa = &sa(t);
    // let lcp = &lcp(t, sa);

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
    let r = r(bwt);
    eprintln!("{name:>40}  {n:>4}  {r:>4}");

    // // slow
    // let (delta, delta_k) = delta(t);
    // let delta_lg = (delta * (n as f32 / delta).log2()) as usize;
    // let delta = delta as usize;
    // let w = w(t, sa, lcp);
    // let chi = chi(t, sa, lcp, print && false);
    // let chi_pd = chi_pd(t, sa, lcp);
    // let chi_pd2 = chi_pd2(t, sa, lcp);
    // // slow

    // let stpd_pos_minus = stpd_pos_minus(t, sa, lcp);
    // // let stpd_pos_plus = stpd_pos_plus(t, sa, lcp);
    // let stpd_lex_minus = stpd_lex_minus(t, sa, lcp);
    // // let stpd_lex_plus = stpd_lex_plus(t, sa, lcp);
    // // let stpd_colex_minus = stpd_colex_minus(t, sa, lcp);
    // // let stpd_colex_plus = stpd_colex_plus(t, sa, lcp);
    // // let stpd_rand = stpd_rand(t, sa, lcp);

    // // let plcp = plcp(t, sa, lcp);

    // // {nodes:>5} {edges:>5}  \
    // // {avg_node_depth:>4} {avg_edge_depth:>4} {inv_avg_node_depth:>5} {inv_avg_edge_depth:>5}  {normalized_tree_size:>4}  \
    // // eprintln!("{name:>40}  {n:>4}  \
    // // {r:>4} {delta:>3} {delta_k:>3} {w:>5} {chi:>5} {chi_pd:>5} {chi_pd2:>5}  {delta_lg:>8} \
    // // {stpd_pos_minus:>5} {stpd_pos_plus:>5} {stpd_lex_minus:>5} {stpd_lex_plus:>5} {stpd_colex_minus:>5} {stpd_colex_plus:>5} {stpd_rand:>5}");
}

fn stpd() {
    env_logger::Builder::from_default_env()
        .format_timestamp_micros()
        .init();
    let t = b"BANANABAANNABBBAAANNANBANANANBANANABANNNAAANABNANNAANABBANNANA";
    // let t = b"AAAAAAAAABAAAAAAA";
    // let t = &[t.as_slice(); 5].concat();
    let t = &relative(1_000_000, 4, 100, 0.001).1;
    // let t = b"BANANABBNNAABBNANAANNABBBAAANNANBANANANBANANABANNNAAANABNANNAANABBANNANAXABBBABABCBANBNANBANANABANAANNANBANABBABANANABNABANNAABBANA";
    // let t = b"AABBCABCBCBB";
    stpd::Stpd::new(t);

    // RopeBWT: 65h for 320 copies of 3.2Gbp => 4.2 Mbp/s  many threads (?)
    // us:      16s for 100 copies of 1  Mbp => 6.2 Mbp/s  1 thread
    // us:     410s for 100 copies of 10 Mbp => 2.4 Mbp/s  1 thread
}

fn main() {
    return stpd();

    // newtest();

    // return;

    // stats(&terminate(fib(6)), true);
    // return;

    header();
    let repeated = relative(200, 4, 20, 0.05);
    let texts = [
        // variants(fib(15)),
        vec![
            // random(3_200_000_00, 4),
            // thue_morse(10),
            // random(1000, 2),
            // relative(500, 2, 2, 0.05),
            // relative(250, 2, 4, 0.05),
            // relative(100, 2, 10, 0.05),
            // relative(50, 2, 20, 0.05),
            // relative(25, 2, 40, 0.05),
            relative(100000, 4, 1, 0.00),
            relative(100000, 4, 2, 0.00),
            relative(100000, 4, 3, 0.00),
            relative(100000, 4, 4, 0.00),
            relative(100000, 4, 1, 0.001),
            relative(100000, 4, 2, 0.001),
            relative(100000, 4, 3, 0.001),
            relative(100000, 4, 4, 0.001),
            relative(100000, 4, 1, 0.01),
            relative(100000, 4, 2, 0.01),
            relative(100000, 4, 3, 0.01),
            relative(100000, 4, 4, 0.01),
            relative(1000000, 4, 1, 0.001),
            relative(1000000, 4, 2, 0.001),
            relative(1000000, 4, 3, 0.001),
            relative(1000000, 4, 4, 0.001),
            relative(1000000, 4, 1, 0.01),
            relative(1000000, 4, 2, 0.01),
            relative(1000000, 4, 3, 0.01),
            relative(1000000, 4, 4, 0.01),
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
        stats(&t, false);
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
