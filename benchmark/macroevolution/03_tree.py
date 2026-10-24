#!/usr/bin/env python3
"""Species on the GTDB r220 bacterial tree: mapping, pruning, parsimony and distance to the references.

Inputs: ``species_counts.npz`` and ``species.tsv`` (02_species_profiles.py), the GTDB r220 bac120 metadata
(species -> representative accession, NCBI taxids) and the GTDB-Tk r220 bac120 reference tree
(``gtdb_r220_bac120_decorated_unrooted.tree``, one tip per species representative), the Swiss-Prot
training taxonomy (``train_taxonomy.tsv``: protein, NCBI taxid; the proteins with experimental GO
annotations that the Stage 0 DIAMOND channel searches), and the GO vocabulary.

Steps
1. Map every species name to its GTDB representative accession; keep the species that are tips.
2. Prune the tree to those tips (branch lengths of collapsed unary nodes are summed), root it at the
   midpoint of its longest tip-to-tip path, write ``tree_pruned.nwk`` and ``tips.tsv``.
3. Reference tips: species that have at least --min-ref-proteins experimentally annotated Swiss-Prot
   proteins (taxid matched to GTDB genomes by ncbi_taxid or ncbi_species_taxid). For every tip the
   patristic distance to the nearest reference tip (``dist_to_ref`` in tips.tsv; 0 for references).
4. Tip states per configuration: a term is present in a species when at least --majority of its
   genomes have it. Fitch parsimony over all terms at once: minimum number of changes per term
   (root-independent), then one most-parsimonious assignment (parent state kept when ambiguous,
   root takes the majority tip state when ambiguous) that gives gains and losses per term and per
   branch. Writes ``parsimony_<cfg>.tsv`` and ``branch_events_<cfg>.tsv`` plus ``nodes.tsv``.
"""
import argparse, csv, gzip, os, sys
from collections import Counter, defaultdict
import numpy as np

sys.setrecursionlimit(1000000)
CFGS = ("noisy_or", "gspa_full")


# ----------------------------------------------------------------------------------------- Newick
def parse_newick(text):
    """Returns (children, parent, length, label) arrays; node 0 is the root of the parse."""
    children, parent, length, label = [[]], [-1], [0.0], [""]
    stack = [0]; cur = 0
    i, n = 0, len(text)
    def new_node(par):
        children.append([]); parent.append(par); length.append(0.0); label.append("")
        children[par].append(len(children) - 1)
        return len(children) - 1
    tokens = []
    buf = []
    quote = False
    for ch in text:
        if ch == "'":
            quote = not quote; buf.append(ch); continue
        if quote:
            buf.append(ch); continue
        if ch in "(),;:":
            if buf:
                tokens.append("".join(buf)); buf = []
            tokens.append(ch)
        elif ch.isspace():
            continue
        else:
            buf.append(ch)
    node = 0
    expect_len = False
    for t in tokens:
        if t == "(":               # the current node is the group; open its first child
            node = new_node(node)
        elif t == ",":
            node = new_node(parent[node])
        elif t == ")":
            node = parent[node] if parent[node] >= 0 else node
        elif t == ":":
            expect_len = True
        elif t == ";":
            break
        else:
            if expect_len:
                length[node] = float(t); expect_len = False
            else:
                label[node] = t.strip("'")
    return children, parent, length, label, 0


def postorder(children, root):
    order, stack = [], [root]
    while stack:
        n = stack.pop(); order.append(n); stack.extend(children[n])
    return order[::-1]


def prune(children, parent, length, label, root, keep):
    """Keep only the tips in `keep`; collapse unary nodes; return new arrays (re-indexed)."""
    order = postorder(children, root)
    alive = {}
    for n in order:
        if not children[n]:
            alive[n] = label[n] in keep
        else:
            alive[n] = any(alive[c] for c in children[n])
    # rebuild
    new_children, new_parent, new_length, new_label = [], [], [], []
    def build(n, par, extra_len):
        kids = [c for c in children[n] if alive[c]]
        if len(kids) == 1 and n != root:
            return build(kids[0], par, extra_len + length[n])
        idx = len(new_children)
        new_children.append([]); new_parent.append(par); new_length.append(length[n] + extra_len if par >= 0 else 0.0)
        new_label.append(label[n])
        for c in kids:
            ci = build(c, idx, 0.0)
            new_children[idx].append(ci)
        return idx
    # if root has a single alive child, descend
    r = root
    while True:
        kids = [c for c in children[r] if alive[c]]
        if len(kids) == 1:
            r = kids[0]
        else:
            break
    build(r, -1, 0.0)
    return new_children, new_parent, new_length, new_label, 0


def to_newick(children, length, label, root):
    out = []
    def rec(n):
        if children[n]:
            out.append("("); first = True
            for c in children[n]:
                if not first: out.append(",")
                first = False; rec(c)
            out.append(")")
        lab = label[n]
        if lab:
            out.append("'" + lab + "'" if any(ch in lab for ch in " ,:;()") else lab)
        if n != root:
            out.append(f":{length[n]:.9g}")
    rec(root)
    return "".join(out) + ";"


def midpoint_root(children, parent, length, label, root):
    """Root at the midpoint of the longest tip-to-tip path (graph two-sweep)."""
    n = len(children)
    adj = [[] for _ in range(n)]
    for v in range(n):
        if parent[v] >= 0:
            adj[v].append((parent[v], length[v])); adj[parent[v]].append((v, length[v]))
    def sweep(src):
        dist = np.full(n, -1.0); prev = np.full(n, -1, dtype=np.int64)
        dist[src] = 0.0; stack = [src]
        while stack:
            u = stack.pop()
            for v, w in adj[u]:
                if dist[v] < 0:
                    dist[v] = dist[u] + w; prev[v] = u; stack.append(v)
        return dist, prev
    tips = [v for v in range(n) if not children[v]]
    d0, _ = sweep(tips[0]); a = int(np.argmax(np.where([not children[v] for v in range(n)], d0, -1)))
    da, prev = sweep(a); b = int(np.argmax(np.where([not children[v] for v in range(n)], da, -1)))
    half = da[b] / 2.0
    # walk from b towards a until passing the midpoint
    v = b
    while da[v] > half:
        u = int(prev[v])
        if da[u] <= half:
            break
        v = u
    u = int(prev[v])  # edge (u, v) contains the midpoint, da[u] <= half < da[v]
    w = da[v] - da[u]; from_u = half - da[u]
    # new root node on edge u-v: re-root the graph from it
    newroot = n
    adj.append([(u, from_u), (v, w - from_u)])
    adj[u] = [(x, l) for x, l in adj[u] if x != v] + [(newroot, from_u)]
    adj[v] = [(x, l) for x, l in adj[v] if x != u] + [(newroot, w - from_u)]
    label = list(label) + [""]
    nc = [[] for _ in range(n + 1)]; npar = [-1] * (n + 1); nlen = [0.0] * (n + 1)
    seen = [False] * (n + 1); seen[newroot] = True; stack = [newroot]
    while stack:
        x = stack.pop()
        for y, l in adj[x]:
            if not seen[y]:
                seen[y] = True; npar[y] = x; nlen[y] = l; nc[x].append(y); stack.append(y)
    # drop any unary node introduced at the old root
    return prune(nc, npar, nlen, label, newroot, {label[v] for v in range(n) if not children[v]})


def nearest_reference(children, parent, length, root, is_ref):
    n = len(children); INF = float("inf")
    down = np.full(n, INF); order = postorder(children, root)
    for v in order:
        if not children[v]:
            down[v] = 0.0 if is_ref[v] else INF
        else:
            down[v] = min(down[c] + length[c] for c in children[v])
    best = down.copy()
    for v in order[::-1]:  # preorder
        p = parent[v]
        if p >= 0:
            # best via parent: parent's best excluding the subtree of v
            alt = best[p] if not np.isclose(best[p], down[v] + length[v]) else INF
            # recompute properly: min over siblings' down + their length, and parent's up value
            sib = [down[c] + length[c] for c in children[p] if c != v]
            up_p = up[p] if p != root else INF
            cand = min(sib + [up_p])
            up[v] = cand + length[v]
            best[v] = min(down[v], up[v])
    return best


def nearest_reference2(children, parent, length, root, is_ref):
    n = len(children); INF = float("inf")
    down = np.full(n, INF); up = np.full(n, INF)
    order = postorder(children, root)
    for v in order:
        if not children[v]:
            down[v] = 0.0 if is_ref[v] else INF
        else:
            down[v] = min(down[c] + length[c] for c in children[v])
    for v in order[::-1]:
        p = parent[v]
        if p < 0:
            continue
        sib = [down[c] + length[c] for c in children[p] if c != v]
        up[v] = min(sib + [up[p]]) + length[v]
    return np.minimum(down, up)


def fitch(children, parent, root, tip_state, tip_of_node, majority_state):
    """tip_state: (n_tips, n_terms) bool. Returns min changes, gains, losses per term, per-branch counts."""
    n = len(children); order = postorder(children, root); n_terms = tip_state.shape[1]
    can0 = np.zeros((n, n_terms), dtype=bool); can1 = np.zeros((n, n_terms), dtype=bool)
    changes = np.zeros(n_terms, dtype=np.int32)
    for v in order:
        if not children[v]:
            s = tip_state[tip_of_node[v]]
            can1[v] = s; can0[v] = ~s
        else:
            i0 = np.ones(n_terms, bool); i1 = np.ones(n_terms, bool)
            u0 = np.zeros(n_terms, bool); u1 = np.zeros(n_terms, bool)
            for c in children[v]:
                i0 &= can0[c]; i1 &= can1[c]; u0 |= can0[c]; u1 |= can1[c]
            empty = ~(i0 | i1)
            can0[v] = np.where(empty, u0, i0); can1[v] = np.where(empty, u1, i1)
            changes += empty
    state = np.zeros((n, n_terms), dtype=bool)
    gains = np.zeros(n_terms, np.int32); losses = np.zeros(n_terms, np.int32)
    br_gain = np.zeros(n, np.int32); br_loss = np.zeros(n, np.int32)
    amb = can0[root] & can1[root]
    state[root] = np.where(amb, majority_state, can1[root])
    for v in order[::-1]:
        if v == root:
            continue
        p = parent[v]; ps = state[p]
        keep = np.where(ps, can1[v], can0[v])
        sv = np.where(keep, ps, can1[v] & ~can0[v])
        state[v] = sv
        g = sv & ~ps; l = ps & ~sv
        gains += g; losses += l; br_gain[v] = g.sum(); br_loss[v] = l.sum()
    return changes, gains, losses, br_gain, br_loss, state


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True, help="dir with species_counts.npz, species.tsv")
    ap.add_argument("--gtdb-metadata", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--train-taxonomy", required=True)
    ap.add_argument("--go-vocab", required=True)
    ap.add_argument("--obo", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--majority", type=float, default=0.5)
    ap.add_argument("--min-ref-proteins", type=int, default=20)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    # ---- GTDB metadata
    sp2rep, taxid2sp, rep_meta = {}, defaultdict(set), {}
    with gzip.open(a.gtdb_metadata, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            tax = r["gtdb_taxonomy"]; sp = tax.split(";")[-1][3:]
            if r["gtdb_representative"] == "t":
                sp2rep[sp] = r["accession"]
                rep_meta[r["accession"]] = (tax, r.get("genome_size", ""), r.get("protein_count", ""))
            for k in ("ncbi_taxid", "ncbi_species_taxid"):
                if r.get(k):
                    taxid2sp[r[k]].add(sp)
    print(f"GTDB r220: {len(sp2rep)} species representatives", flush=True)
    # ---- our species
    d = np.load(os.path.join(a.profiles, "species_counts.npz"))
    species = list(d["species"]); n_gen = d["n_genomes"]
    sp_rows = {r["species"]: r for r in csv.DictReader(open(os.path.join(a.profiles, "species.tsv")), delimiter="\t")}
    mapped = {sp: sp2rep[sp] for sp in species if sp in sp2rep}
    unmapped = [sp for sp in species if sp not in sp2rep]
    print(f"{len(mapped)} of {len(species)} species map to a GTDB r220 representative; unmapped genomes: "
          f"{int(sum(n_gen[i] for i, sp in enumerate(species) if sp not in sp2rep))}", flush=True)
    with open(os.path.join(a.out, "unmapped_species.tsv"), "w") as fh:
        fh.write("species\tn_genomes\n")
        for i, sp in enumerate(species):
            if sp not in sp2rep:
                fh.write(f"{sp}\t{n_gen[i]}\n")
    # ---- tree
    text = open(a.tree).read()
    ch, par, ln, lab, root = parse_newick(text)
    tips_all = {lab[v] for v in range(len(ch)) if not ch[v]}
    print(f"tree: {len(ch)} nodes, {len(tips_all)} tips", flush=True)
    keep = set(mapped.values()) & tips_all
    missing = set(mapped.values()) - tips_all
    print(f"{len(keep)} representatives are tips; {len(missing)} not in the tree", flush=True)
    ch, par, ln, lab, root = prune(ch, par, ln, lab, root, keep)
    ch, par, ln, lab, root = midpoint_root(ch, par, ln, lab, root)
    n = len(ch)
    tips = [v for v in range(n) if not ch[v]]
    print(f"pruned, midpoint-rooted: {n} nodes, {len(tips)} tips", flush=True)
    open(os.path.join(a.out, "tree_pruned.nwk"), "w").write(to_newick(ch, ln, lab, root) + "\n")
    rep2sp = {v: k for k, v in mapped.items()}
    tip_sp = [rep2sp[lab[v]] for v in tips]
    # ---- reference tips (Swiss-Prot experimental proteins)
    prot_per_taxid = Counter()
    with open(a.train_taxonomy) as fh:
        for line in fh:
            p, t = line.rstrip("\n").split("\t")[:2]
            prot_per_taxid[t] += 1
    ref_prot = Counter()
    for t, c in prot_per_taxid.items():
        for sp in taxid2sp.get(t, ()):
            ref_prot[sp] += c
    is_ref = np.zeros(n, bool)
    for v, sp in zip(tips, tip_sp):
        is_ref[v] = ref_prot.get(sp, 0) >= a.min_ref_proteins
    print(f"{int(is_ref.sum())} reference tips (>= {a.min_ref_proteins} experimental Swiss-Prot proteins)", flush=True)
    dist = nearest_reference2(ch, par, ln, root, is_ref)
    depth = np.zeros(n)
    for v in postorder(ch, root)[::-1]:
        if par[v] >= 0: depth[v] = depth[par[v]] + ln[v]
    # ---- tip states
    with open(a.go_vocab) as fh:
        vocab = [l.rstrip("\n") for l in fh]
    names = {}
    cur = None
    for line in open(a.obo):
        if line.startswith("id: GO:"): cur = line[4:].strip()
        elif line.startswith("name: ") and cur: names[cur] = line[6:].strip(); cur = None
    sidx = {sp: i for i, sp in enumerate(species)}
    tip_rows = np.array([sidx[sp] for sp in tip_sp])
    tip_of_node = np.full(n, -1); tip_of_node[tips] = np.arange(len(tips))
    with open(os.path.join(a.out, "tips.tsv"), "w") as fh:
        fh.write("node\tspecies\taccession\tgtdb_taxonomy\tn_genomes\tis_reference\tn_ref_proteins\tdist_to_ref\tdepth\n")
        for v, sp in zip(tips, tip_sp):
            fh.write(f"{v}\t{sp}\t{lab[v]}\t{rep_meta[lab[v]][0]}\t{n_gen[sidx[sp]]}\t{int(is_ref[v])}\t{ref_prot.get(sp, 0)}\t{dist[v]:.6f}\t{depth[v]:.6f}\n")
    with open(os.path.join(a.out, "nodes.tsv"), "w") as fh:
        fh.write("node\tparent\tlength\tlabel\tn_tips\tdepth\n")
        ntips = np.zeros(n, np.int64)
        for v in postorder(ch, root):
            ntips[v] = 1 if not ch[v] else sum(ntips[c] for c in ch[v])
        for v in range(n):
            fh.write(f"{v}\t{par[v]}\t{ln[v]:.9g}\t{lab[v]}\t{ntips[v]}\t{depth[v]:.6f}\n")
    for cfg in CFGS:
        counts = d[f"{cfg}_counts"]
        frac = counts[tip_rows].astype(np.float32) / n_gen[tip_rows][:, None]
        tip_state = frac >= a.majority
        n_present = tip_state.sum(0)
        informative = (n_present >= 1) & (n_present <= len(tips) - 1)
        terms_idx = np.where(informative)[0]
        ts = tip_state[:, terms_idx]
        majority_state = n_present[terms_idx] * 2 > len(tips)
        print(f"[{cfg}] {len(terms_idx)} variable terms over {len(tips)} tips; running Fitch", flush=True)
        changes, gains, losses, br_gain, br_loss, state = fitch(ch, par, root, ts, tip_of_node, majority_state)
        with open(os.path.join(a.out, f"parsimony_{cfg}.tsv"), "w") as fh:
            fh.write("go_id\tname\tn_tips_present\tmin_changes\tgains\tlosses\tconsistency_index\troot_state\n")
            for j, t in enumerate(terms_idx):
                ci = 1.0 / changes[j] if changes[j] else 1.0
                fh.write(f"{vocab[t]}\t{names.get(vocab[t], '')}\t{int(n_present[t])}\t{int(changes[j])}\t{int(gains[j])}\t{int(losses[j])}\t{ci:.4f}\t{int(state[root, j])}\n")
        with open(os.path.join(a.out, f"branch_events_{cfg}.tsv"), "w") as fh:
            fh.write("node\tparent\tlength\tn_tips\tn_gains\tn_losses\tn_present\n")
            for v in range(n):
                fh.write(f"{v}\t{par[v]}\t{ln[v]:.9g}\t{ntips[v]}\t{br_gain[v]}\t{br_loss[v]}\t{int(state[v].sum())}\n")
        # ancestral states of the root and of every internal node with a taxon label
        np.savez_compressed(os.path.join(a.out, f"ancestral_{cfg}.npz"), terms=np.array([vocab[t] for t in terms_idx]),
                            state=np.packbits(state, axis=1), n_terms=np.int64(len(terms_idx)))
        print(f"[{cfg}] total min changes {int(changes.sum())}, gains {int(gains.sum())}, losses {int(losses.sum())}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
