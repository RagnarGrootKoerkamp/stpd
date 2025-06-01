use chi::*;

fn header() {
    let name = "name";
    let n = "n";
    let r = "r";
    let delta = "δ";
    let delta_k = "δₖ";
    let w = "W";
    let chi = "χ";
    let chi_pd = "χ pd";
    let chi_pd2 = "χ set-pd";
    let nodes = "nodes";
    let edges = "edges";
    let avg_node_depth = "n.dp";
    let avg_edge_depth = "e.dp";
    let inv_avg_node_depth = "in.dp";
    let inv_avg_edge_depth = "ie.dp";
    let normalized_tree_size = "ntsz";
    eprintln!("{name:>30}  {n:>4}  {nodes:>5} {edges:>5}  {avg_node_depth:>4} {avg_edge_depth:>4} {inv_avg_node_depth:>5} {inv_avg_edge_depth:>5}  {normalized_tree_size:>4}  {r:>4} {delta:>3} {delta_k:>3} {w:>5} {chi:>5} {chi_pd:>5} {chi_pd2:>5}");
}

fn stats((name, t): &(String, T), print: bool) {
    let n = t.len();
    let sa = &sa(t);
    let lcp = &lcp(t, sa);

    let mut nodes = 0;
    let mut node_depth = 0;
    for n in tree_nodes(t, sa, lcp) {
        nodes += 1;
        node_depth += n.len();
    }
    let avg_node_depth = node_depth / nodes;
    let inv_avg_node_depth = n * nodes / node_depth;

    let mut edges = 0;
    let mut edge_depth = 0;
    for e in tree_edges(t, sa, lcp) {
        edges += 1;
        edge_depth += e.len();
    }
    let avg_edge_depth = edge_depth / edges;
    let inv_avg_edge_depth = n * edges / edge_depth;

    let mut tree_size = n * (n - 1) / 2;
    for n in tree_nodes(t, sa, lcp) {
        tree_size += n.len() + 1;
    }
    for e in tree_edges(t, sa, lcp) {
        tree_size -= e.len();
    }
    let normalized_tree_size = tree_size / n;

    let bwt = &bwt(t, sa);
    let r = r(bwt);
    let (delta, delta_k) = delta(t);
    let w = w(t, sa, lcp);
    let chi = chi(t, sa, lcp, print);
    let chi_pd = chi_pd(t, sa, lcp);
    let chi_pd2 = chi_pd2(t, sa, lcp);

    eprintln!("{name:>30}  {n:>4}  {nodes:>5} {edges:>5}  {avg_node_depth:>4} {avg_edge_depth:>4} {inv_avg_node_depth:>5} {inv_avg_edge_depth:>5}  {normalized_tree_size:>4}  {r:>4} {delta:>3.0} {delta_k:>3} {w:>5} {chi:>5} {chi_pd:>5} {chi_pd2:>5}");
}

fn main() {
    header();
    use chi::strings::*;
    let texts = [
        // variants(fib(6)),
        variants(fib(15)),
        vec![thue_morse(10)],
        vec![random(1000, 2)],
        vec![random(1000, 4)],
        vec![relative(100, 2, 10, 0.05)],
        vec![relative(100, 4, 10, 0.05)],
    ]
    .concat();
    for t in texts {
        stats(&t, false);
    }

    // stats(&random(1000, 2), true);
}
