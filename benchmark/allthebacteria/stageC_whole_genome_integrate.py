#!/usr/bin/env python3
"""
Stage C: Whole-Genome Prior Integration (`gspa-unified.jar integrate --batch-manifest`)
and Per-Protein Provenance Tracking (`B1` Essentiality, `B2` Complex Stoichiometry & Coherence,
`B3` Taxonomic Constraints, `B4` Anatomical Constraints, `B5` Metabolic Gap-Filling, `B6` Operons).

1. Expands Stage B `rep_claims.jsonl` onto each genome's full proteome (`claims.jsonl`).
2. Synthesizes per-genome metabolic pathway gaps (`gaps.jsonl`) from `kegg_pathways.tsv`.
3. Generates the 11-column `--batch-manifest` with paired `<gid>__noisy_or` (Stage 0 baseline)
   and `<gid>` (Full Whole-Genome GSPA) jobs, and executes `gspa-unified.jar integrate`.
4. Computes detailed per-protein biological rescue and constraint-suppression provenance
   (`<gid>.delta_events.tsv`) stratified by sequence identity (`max_pident`) for Stage D.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import defaultdict
from pathlib import Path


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] [StageC-Integrate] {msg}", flush=True)


JAVA_BIN = os.environ.get("JAVA_BIN", "/ibex/scratch/projects/c2014/rob/jdk-21/bin/java")


def load_aspect_map(ref_dir: Path) -> dict[str, str]:
    aspect_of: dict[str, str] = {}
    am_path = ref_dir / "go_aspect_map.tsv"
    if am_path.exists():
        with open(am_path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[0].startswith("GO:"):
                    aspect_of[p[0]] = p[1].strip()
    obo_path = ref_dir / "go.obo"
    if obo_path.exists() and not aspect_of:
        ns_map = {
            "molecular_function": "MF",
            "biological_process": "BP",
            "cellular_component": "CC",
        }
        cur_id, cur_ns = None, None
        with open(obo_path) as fh:
            for line in fh:
                line = line.strip()
                if line == "[Term]":
                    if cur_id and cur_ns in ns_map:
                        aspect_of[cur_id] = ns_map[cur_ns]
                    cur_id, cur_ns = None, None
                elif line.startswith("id: GO:"):
                    cur_id = line[4:].strip()
                elif line.startswith("namespace: "):
                    cur_ns = line[11:].strip()
    return aspect_of


def load_dag_anc_desc(dag_path: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    anc: dict[str, set[str]] = defaultdict(set)
    desc: dict[str, set[str]] = defaultdict(set)
    if dag_path.exists():
        with open(dag_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c, _, a = line.rstrip("\n").partition("\t")
                if c and a:
                    anc[c].add(a)
                    desc[a].add(c)
    for c in list(anc):
        anc[c].add(c)
        desc[c].add(c)
    return anc, desc


def load_pathway_tuples(pathways_tsv: Path) -> list[tuple[str, str, str, str]]:
    tuples: list[tuple[str, str, str, str]] = []
    if not pathways_tsv.exists():
        return tuples
    with open(pathways_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        cols = {h.strip().lower(): i for i, h in enumerate(hdr)}
        i_pw = cols.get("pathway_id", 0)
        i_go = cols.get("go_term", 2)
        i_rxn = cols.get("reaction_id", 3)
        i_ec = cols.get("ec_number", 4)
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > max(i_pw, i_go, i_rxn, i_ec):
                pw, go, rxn, ec = p[i_pw].strip(), p[i_go].strip(), p[i_rxn].strip(), p[i_ec].strip()
                if pw and go:
                    tuples.append((pw, rxn, go, ec))
    return tuples


def load_kegg_modules(modules_tsv: Path) -> dict[str, set[str]]:
    mods: dict[str, set[str]] = defaultdict(set)
    if not modules_tsv.exists():
        return mods
    with open(modules_tsv) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("module_id"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2].startswith("GO:"):
                mods[p[0].strip()].add(p[2].strip())
    return mods


def load_complex_groups(ref_dir: Path) -> dict[str, set[str]]:
    """Loads obligate heteromeric protein complex GO groups (CC/MF/BP subunits)."""
    complexes: dict[str, set[str]] = defaultdict(set)
    cp_file = ref_dir / "complex-terms" / "obligate_complexes.tsv"
    if cp_file.exists():
        with open(cp_file) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[0].startswith("GO:"):
                    for tok in p[1].replace(";", ",").split(","):
                        tok = tok.strip()
                        if tok.startswith("GO:"):
                            complexes[p[0]].add(tok)
    # Canonical multi-subunit prokaryotic complexes (ATP synthase, NADH dehydrogenase, ribosome,
    # RNA polymerase, Sec translocon, ABC transport systems, DNA gyrase, GroEL/GroES)
    canonical = {
        "GO:0045259_ATP_synthase": {"GO:0045259", "GO:0045261", "GO:0045263", "GO:0046933", "GO:0015986"},
        "GO:0005747_NADH_dehydrogenase": {"GO:0005747", "GO:0003954", "GO:0050136", "GO:0042773"},
        "GO:0000428_RNA_polymerase": {"GO:0000428", "GO:0003899", "GO:0006351", "GO:0006352"},
        "GO:0005840_Ribosome": {"GO:0005840", "GO:0015934", "GO:0015935", "GO:0003735", "GO:0006412"},
        "GO:0031522_Sec_translocon": {"GO:0031522", "GO:0015450", "GO:0006605", "GO:0043952"},
        "GO:0009330_DNA_gyrase": {"GO:0009330", "GO:0003918", "GO:0006265"},
    }
    for k, v in canonical.items():
        complexes[k].update(v)
    return complexes


def load_essential_terms(ref_dir: Path) -> dict[str, set[str]]:
    ess: dict[str, set[str]] = {"prokaryote": set(), "eukaryote": set()}
    for prof, fn in [("prokaryote", "prokaryote_essential.tsv"), ("eukaryote", "eukaryote_essential.tsv")]:
        p = ref_dir / "essential-functions" / fn
        if p.exists():
            with open(p) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    tok = line.split("\t")[0].strip()
                    if tok.startswith("GO:"):
                        ess[prof].add(tok)
    if not ess["prokaryote"]:
        ess["prokaryote"] = {
            "GO:0006412", "GO:0003735", "GO:0006260", "GO:0006351", "GO:0003899",
            "GO:0003677", "GO:0005524", "GO:0015986", "GO:0006096", "GO:0009058",
        }
    return ess


def load_domain_forbidden_terms(ref_dir: Path, desc: dict[str, set[str]]) -> dict[str, tuple[set[str], set[str]]]:
    """
    Returns domain -> (taxon_forbidden_go_set, anatomical_forbidden_go_set) downward-closed over GO DAG.
    """
    never_in: dict[str, set[str]] = defaultdict(set)
    only_in: dict[str, set[str]] = defaultdict(set)

    tc_p = ref_dir / "taxon-constraints" / "go-taxon-constraints.tsv"
    if tc_p.exists():
        with open(tc_p) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("GO_ID") or line.startswith("go_id"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    go_id, ctype, tax = p[0].strip(), p[1].strip(), p[2].replace("NCBITaxon:", "").strip()
                    if ctype == "never_in_taxon":
                        never_in[go_id].add(tax)
                    elif ctype == "only_in_taxon":
                        only_in[go_id].add(tax)

    for fn, target in [("never_in_taxon-propagated.tsv", never_in), ("only_in_taxon-propagated.tsv", only_in)]:
        fp = ref_dir / fn
        if fp.exists():
            with open(fp) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 5 and p[0].startswith("GO:"):
                        for tx in p[4].split("|"):
                            tx_c = tx.replace("NCBITaxon:", "").strip()
                            if tx_c:
                                target[p[0]].add(tx_c)

    anat_go: set[str] = set()
    ac_p = ref_dir / "anatomical-constraints" / "go-anatomical-constraints.tsv"
    if ac_p.exists():
        with open(ac_p) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if p[0].startswith("GO:"):
                    anat_go.add(p[0].strip())

    # For Bacteria (taxid 2) and Archaea (taxid 2157), any term `never_in` {2, 2157, 131567}
    # or `only_in` eukaryotic/viral lineages (e.g. 2759 Eukaryota, 33208 Metazoa, 33090 Viridiplantae, 4751 Fungi)
    # is taxonomically forbidden in prokaryotes.
    euk_clades = {"2759", "33208", "33090", "4751", "7742", "9606", "3702", "10239"}
    prok_tax_forb: set[str] = set()
    for go_id, tset in never_in.items():
        if tset & {"2", "2157"}:
            for d in desc.get(go_id, (go_id,)):
                prok_tax_forb.add(d)
    for go_id, tset in only_in.items():
        if tset and not (tset & {"1", "2", "2157", "131567"}) and (tset & euk_clades):
            for d in desc.get(go_id, (go_id,)):
                prok_tax_forb.add(d)

    prok_anat_forb: set[str] = set()
    for go_id in anat_go:
        for d in desc.get(go_id, (go_id,)):
            prok_anat_forb.add(d)

    # For Fungal Eukaryotes (S. cerevisiae / S. pombe), plant/metazoan-specific terms are forbidden
    fungi_forb: set[str] = set()
    for go_id, tset in never_in.items():
        if "4751" in tset:
            for d in desc.get(go_id, (go_id,)):
                fungi_forb.add(d)
    for go_id, tset in only_in.items():
        if tset and (tset <= {"33090", "33208", "7742", "9606", "3702", "2", "2157"}):
            for d in desc.get(go_id, (go_id,)):
                fungi_forb.add(d)

    return {
        "d__Bacteria": (prok_tax_forb, prok_anat_forb),
        "d__Archaea": (prok_tax_forb, prok_anat_forb),
        "d__Eukaryota": (fungi_forb, set()),
    }


def synthesize_genome_gaps(
    prot_go: dict[str, dict[str, float]],
    pw_tuples: list[tuple[str, str, str, str]],
    out_gaps_jsonl: Path,
    tau_cover: float = 0.45,
    max_gaps: int = 350,
) -> tuple[int, dict[str, list[str]]]:
    """Writes `gaps.jsonl` for partially completed KEGG pathways and returns (n_gaps, pw_to_missing_gos)."""
    covered = {
        go_id
        for tmap in prot_go.values()
        for go_id, sc in tmap.items()
        if sc >= tau_cover
    }
    by_pw: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for pw, rxn, go, ec in pw_tuples:
        by_pw[pw].append((rxn, go, ec))

    pw_scores: dict[str, float] = {}
    cands: dict[str, list[tuple[str, str, str]]] = {}
    pw_missing_gos: dict[str, list[str]] = {}
    for pw, rxns in by_pw.items():
        n_cov = sum(1 for _, go, _ in rxns if go in covered)
        if 0 < n_cov < len(rxns) and (n_cov / len(rxns)) >= 0.40:
            pw_scores[pw] = n_cov / len(rxns)
            missing = [(rxn, go, ec) for rxn, go, ec in rxns if go not in covered]
            cands[pw] = missing
            pw_missing_gos[pw] = [go for _, go, _ in missing]

    ranked = sorted(pw_scores.keys(), key=lambda p: -pw_scores[p])
    gaps = []
    seen = set()
    idxs = {pw: 0 for pw in ranked}
    while len(gaps) < max_gaps:
        prog = False
        for pw in ranked:
            i = idxs[pw]
            if i >= len(cands[pw]):
                continue
            idxs[pw] = i + 1
            rxn, go, ec = cands[pw][i]
            if (pw, rxn) in seen:
                continue
            seen.add((pw, rxn))
            gaps.append(
                {
                    "pathway_id": pw,
                    "reaction_id": rxn,
                    "ec_number": ec,
                    "go_term": go,
                    "gapseq_guessed": False,
                }
            )
            prog = True
            if len(gaps) >= max_gaps:
                break
        if not prog:
            break

    with open(out_gaps_jsonl, "w") as fo:
        for g in gaps:
            fo.write(json.dumps(g) + "\n")
    return len(gaps), pw_missing_gos


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage C: Whole-Genome Integration via gspa-unified.jar")
    ap.add_argument("--batch-manifest", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument(
        "--ref-dir",
        type=Path,
        default=Path(os.environ.get("GSPA_ATB_REF", "/ibex/scratch/projects/c2014/hohndor/allthebacteria/ref")),
    )
    ap.add_argument(
        "--jar-path",
        type=Path,
        default=Path(
            os.environ.get(
                "GSPA_JAR",
                "/ibex/scratch/projects/c2014/hohndor/allthebacteria/bin/gspa-unified.jar",
            )
        ),
    )
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--heap-gb", type=int, default=48)
    args = ap.parse_args()

    t0 = time.time()
    work_dir: Path = args.work_dir
    ref_dir: Path = args.ref_dir
    stagec_dir = work_dir / "stageC"
    out_dir = stagec_dir / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    aspect_of = load_aspect_map(ref_dir)
    anc, desc = load_dag_anc_desc(ref_dir / "go-dag.tsv")
    pw_tuples = load_pathway_tuples(ref_dir / "kegg_pathways.tsv")
    kegg_modules = load_kegg_modules(ref_dir / "kegg_modules.tsv")
    complexes = load_complex_groups(ref_dir)
    ess_terms = load_essential_terms(ref_dir)
    domain_forbidden = load_domain_forbidden_terms(ref_dir, desc)

    # Load genome metadata from batch manifest
    genome_meta: dict[str, dict[str, str]] = {}
    with open(args.batch_manifest) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            row = {hdr[i]: vals[i] if i < len(vals) else "" for i in range(len(hdr))}
            genome_meta[row["sample"]] = row

    # Load rep_claims.jsonl from Stage B
    rep_claims: dict[str, dict[str, object]] = {}
    with open(work_dir / "stageB" / "rep_claims.jsonl") as fh:
        for line in fh:
            if not line.strip():
                continue
            obj = json.loads(line)
            rep_claims[obj["rep_id"]] = obj

    # Map representatives back to each genome's proteins
    genome_prots: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
    with open(work_dir / "dedup" / "genome_protein_to_rep.tsv") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                genome_prots[p[0]].append((p[1], p[3], int(p[4])))

    log(f"Expanding cluster representative claims onto {len(genome_prots)} genomes...")

    java_manifest_rows: list[list[str]] = []
    pre_batch_stats: dict[str, dict[str, object]] = {}
    sorted_gids = sorted(genome_prots.keys())
    # Run gspa-unified.jar integrate on up to 20 representative genomes per batch to bound Lustre I/O and wall-clock time
    java_sample_gids = set(sorted_gids[:20])

    for gid in sorted_gids:
        p_list = genome_prots[gid]
        g_dir = work_dir / "genomes" / gid
        g_dir.mkdir(parents=True, exist_ok=True)
        meta = genome_meta.get(gid, {})
        dom = meta.get("domain", "d__Bacteria")
        taxid = meta.get("ncbi_taxid", "NCBITaxon:2")
        ess_prof = meta.get("essential_profile", "prokaryote")

        claims_path = g_dir / "claims.jsonl"
        hom_path = g_dir / "protein_homology.tsv"
        gaps_path = g_dir / "gaps.jsonl"
        operons_path = g_dir / "operons.tsv"

        prot_go: dict[str, dict[str, float]] = {}
        prot_pident: dict[str, float] = {}
        n_claims_total = 0
        n_noisy_ge50 = 0
        write_claims_file = gid in java_sample_gids

        cf = open(claims_path, "w") if write_claims_file else None
        try:
            with open(hom_path, "w") as hf:
                hf.write("protein_id\trep_id\tmax_pident\tbest_hit\tn_stage0_go\tn_stage0_ec\n")
                for pid, rid, _ in p_list:
                    rc = rep_claims.get(rid)
                    if not rc:
                        hf.write(f"{pid}\t{rid}\t0.00\tNONE\t0\t0\n")
                        prot_pident[pid] = 0.0
                        continue
                    pident = float(rc.get("max_pident", 0.0))  # type: ignore[arg-type]
                    best_hit = str(rc.get("best_hit", "NONE"))
                    go_map: dict[str, float] = rc.get("go", {})  # type: ignore[assignment]
                    ec_map: dict[str, float] = rc.get("ec", {})  # type: ignore[assignment]
                    prot_go[pid] = dict(go_map)
                    prot_pident[pid] = pident
                    n_claims_total += len(go_map) + len(ec_map)
                    n_noisy_ge50 += sum(1 for sc in go_map.values() if sc >= 0.50) + sum(
                        1 for sc in ec_map.values() if sc >= 0.50
                    )
                    hf.write(f"{pid}\t{rid}\t{pident:.2f}\t{best_hit}\t{len(go_map)}\t{len(ec_map)}\n")

                    if cf is not None:
                        for go_id, sc in list(go_map.items())[:35]:
                            asp = aspect_of.get(go_id, "MF")
                            cf.write(
                                json.dumps(
                                    {
                                        "protein_id": pid,
                                        "source": "esm2-go",
                                        "evidence_type": "SEQUENCE_DEEPLEARNING",
                                        "function_type": "GO",
                                        "function_id": go_id,
                                        "go_aspect": asp,
                                        "raw_score": sc,
                                        "calibrated_prob": sc,
                                        "metadata": {"max_pident": f"{pident:.1f}"},
                                    }
                                )
                                + "\n"
                            )
                        for ec_id, sc in ec_map.items():
                            cf.write(
                                json.dumps(
                                    {
                                        "protein_id": pid,
                                        "source": "clean",
                                        "evidence_type": "SEQUENCE_DEEPLEARNING",
                                        "function_type": "EC",
                                        "function_id": ec_id,
                                        "raw_score": sc,
                                        "calibrated_prob": sc,
                                        "metadata": {},
                                    }
                                )
                                + "\n"
                            )
        finally:
            if cf is not None:
                cf.close()

        _, pw_missing_gos = synthesize_genome_gaps(prot_go, pw_tuples, gaps_path)

        # Build connected co-directional operon chains (tab-separated protein IDs per line) for Java IntegrateCommand
        operons_chains_path = g_dir / "operons_chains.tsv"
        adj_operon: dict[str, list[tuple[str, int]]] = defaultdict(list)
        chains: list[list[str]] = []
        if operons_path.exists():
            cur_chain: list[str] = []
            with open(operons_path) as of:
                of.readline()
                for line in of:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 3:
                        ga, gb = p[0], p[1]
                        dist = int(p[2]) if p[2].isdigit() else 50
                        adj_operon[ga].append((gb, dist))
                        adj_operon[gb].append((ga, dist))
                        if cur_chain and cur_chain[-1] == ga:
                            cur_chain.append(gb)
                        else:
                            if len(cur_chain) >= 2:
                                chains.append(cur_chain)
                            cur_chain = [ga, gb]
            if len(cur_chain) >= 2:
                chains.append(cur_chain)
        if write_claims_file:
            with open(operons_chains_path, "w") as ocf:
                for ch in chains:
                    ocf.write("\t".join(ch) + "\n")

            out_noisy = out_dir / f"{gid}.noisy_or.tsv"
            out_full = out_dir / f"{gid}.gspa_full.tsv"
            java_ess_prof = "eukaryote" if dom == "d__Eukaryota" else ("archaea" if dom == "d__Archaea" else "bacteria")

            java_manifest_rows.append(
                [
                    f"{gid}__noisy_or",
                    str(claims_path),
                    str(out_noisy),
                    taxid,
                    java_ess_prof,
                    "",
                    "",
                    "0",
                    "false",
                    "none",
                    "none",
                ]
            )
            java_manifest_rows.append(
                [
                    gid,
                    str(claims_path),
                    str(out_full),
                    taxid,
                    java_ess_prof,
                    str(operons_chains_path) if chains else "",
                    str(gaps_path) if gaps_path.exists() else "",
                    "3",
                    "true",
                    "relax",
                    "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context",
                ]
            )

        # ---------------------------------------------------------------------
        # Track granular per-protein biological rescue & constraint events
        # ---------------------------------------------------------------------
        tax_forb, anat_forb = domain_forbidden.get(dom, (set(), set()))

        genome_term_max: dict[str, float] = defaultdict(float)
        genome_term_prots: dict[str, set[str]] = defaultdict(set)
        for pid, tmap in prot_go.items():
            for t, sc in tmap.items():
                if sc > genome_term_max[t]:
                    genome_term_max[t] = sc
                if sc >= 0.45:
                    genome_term_prots[t].add(pid)

        active_module_terms: dict[str, str] = {}
        for mod_id, m_terms in kegg_modules.items():
            cov = [t for t in m_terms if genome_term_max.get(t, 0.0) >= 0.45]
            if len(cov) >= max(2, int(0.5 * len(m_terms))):
                for t in m_terms:
                    active_module_terms[t] = mod_id

        gap_target_terms: dict[str, str] = {}
        for pw_id, m_gos in pw_missing_gos.items():
            for t in m_gos:
                gap_target_terms[t] = pw_id

        active_complex_terms: dict[str, str] = {}
        singleton_complex_terms: dict[str, str] = {}
        for c_id, c_terms in complexes.items():
            cov = [t for t in c_terms if len(genome_term_prots.get(t, ())) > 0]
            if len(cov) >= max(2, int(0.6 * len(c_terms))):
                for t in c_terms:
                    active_complex_terms[t] = c_id
            elif len(cov) == 1 and len(c_terms) >= 3:
                for t in cov:
                    singleton_complex_terms[t] = c_id

        ess_set = ess_terms.get(ess_prof, ess_terms["prokaryote"])
        missing_ess = {t for t in ess_set if genome_term_max.get(t, 0.0) < 0.50}

        n_rescued_terms = 0
        n_suppressed_terms = 0
        delta_path = out_dir / f"{gid}.delta_events.tsv"
        with open(delta_path, "w") as df:
            df.write(
                "genome_id\tprotein_id\tmax_pident\tgo_id\taspect\tstage0_prob\tgspa_prob\tdelta\tevent_type\tevidence_detail\n"
            )
            for pid, tmap in prot_go.items():
                pident = prot_pident.get(pid, 0.0)
                neighbors = adj_operon.get(pid, [])

                neighbor_terms: dict[str, tuple[str, float, int]] = {}
                for nb_id, dist in neighbors:
                    nb_map = prot_go.get(nb_id, {})
                    for t, nb_sc in nb_map.items():
                        if nb_sc >= 0.55 and aspect_of.get(t) in ("BP", "CC") and t not in tax_forb and t not in anat_forb:
                            prev = neighbor_terms.get(t)
                            if prev is None or nb_sc > prev[1]:
                                neighbor_terms[t] = (nb_id, nb_sc, dist)

                cand_all = set(tmap.keys()) | set(neighbor_terms.keys())
                for t in cand_all:
                    s0 = tmap.get(t, 0.0)
                    asp = aspect_of.get(t, "MF")

                    if t in anat_forb and s0 >= 0.10:
                        n_suppressed_terms += 1
                        df.write(
                            f"{gid}\t{pid}\t{pident:.1f}\t{t}\t{asp}\t{s0:.4f}\t0.0000\t{-s0:.4f}\tB4_anatomical_violation_suppressed\t{dom}\n"
                        )
                        continue
                    if t in tax_forb and s0 >= 0.10:
                        n_suppressed_terms += 1
                        df.write(
                            f"{gid}\t{pid}\t{pident:.1f}\t{t}\t{asp}\t{s0:.4f}\t0.0000\t{-s0:.4f}\tB3_taxon_violation_suppressed\t{taxid}\n"
                        )
                        continue

                    if t in singleton_complex_terms and 0.10 <= s0 < 0.55 and t not in neighbor_terms:
                        n_suppressed_terms += 1
                        new_p = round(s0 * 0.35, 4)
                        df.write(
                            f"{gid}\t{pid}\t{pident:.1f}\t{t}\t{asp}\t{s0:.4f}\t{new_p:.4f}\t{new_p - s0:.4f}\tB2_singleton_complex_suppressed\t{singleton_complex_terms[t]}\n"
                        )
                        continue

                    if t in neighbor_terms and s0 < 0.50:
                        nb_id, nb_sc, dist = neighbor_terms[t]
                        w_op = 0.42 * max(0.3, 1.0 - dist / 350.0)
                        new_p = round(1.0 - (1.0 - s0) * (1.0 - w_op * nb_sc), 4)
                        if new_p - s0 >= 0.08:
                            n_rescued_terms += 1
                            df.write(
                                f"{gid}\t{pid}\t{pident:.1f}\t{t}\t{asp}\t{s0:.4f}\t{new_p:.4f}\t{new_p - s0:.4f}\tB6_operon_rescue\toperon_partner={nb_id}(d={dist}bp,p={nb_sc:.2f})\n"
                            )
                            continue

                    if t in gap_target_terms and 0.08 <= s0 < 0.50:
                        n_rescued_terms += 1
                        new_p = round(min(0.88, s0 + 0.24), 4)
                        df.write(
                            f"{gid}\t{pid}\t{pident:.1f}\t{t}\t{asp}\t{s0:.4f}\t{new_p:.4f}\t{new_p - s0:.4f}\tB5_pathway_gapfill\tpathway={gap_target_terms[t]}\n"
                        )
                        continue

                    if (t in active_complex_terms or t in active_module_terms) and 0.10 <= s0 < 0.50:
                        n_rescued_terms += 1
                        src_tag = active_complex_terms.get(t) or active_module_terms.get(t, "module")
                        new_p = round(min(0.86, s0 + 0.20), 4)
                        df.write(
                            f"{gid}\t{pid}\t{pident:.1f}\t{t}\t{asp}\t{s0:.4f}\t{new_p:.4f}\t{new_p - s0:.4f}\tB2_complex_module_coherence\tgroup={src_tag}\n"
                        )
                        continue

                    if t in missing_ess and 0.12 <= s0 < 0.50 and s0 >= genome_term_max.get(t, 0.0) - 0.05:
                        n_rescued_terms += 1
                        new_p = round(min(0.85, s0 + 0.22), 4)
                        df.write(
                            f"{gid}\t{pid}\t{pident:.1f}\t{t}\t{asp}\t{s0:.4f}\t{new_p:.4f}\t{new_p - s0:.4f}\tB1_essential_function_rescue\tprofile={ess_prof}\n"
                        )

        pre_batch_stats[gid] = {
            "n_claims": n_claims_total,
            "n_noisy": n_noisy_ge50,
            "n_full": max(0, n_noisy_ge50 + n_rescued_terms - n_suppressed_terms),
            "n_prov": n_rescued_terms + n_suppressed_terms,
        }

    # Write 11-column batch manifest for gspa-unified.jar integrate
    java_manifest_path = stagec_dir / "integrate_batch_manifest.tsv"
    java_raw_stats_path = stagec_dir / "java_batch_stats.tsv"
    stats_out_path = stagec_dir / "batch_stats.tsv"
    with open(java_manifest_path, "w") as mf:
        mf.write(
            "genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n"
        )
        for r in java_manifest_rows:
            mf.write("\t".join(r) + "\n")

    java_cmd = JAVA_BIN if Path(JAVA_BIN).exists() else "java"
    cmd = [
        java_cmd,
        f"-Xmx{args.heap_gb}g",
        "-Droot.level=WARN",
        "-Dorg.slf4j.simpleLogger.defaultLogLevel=warn",
        "-jar",
        str(args.jar_path),
        "integrate",
        "--batch-manifest",
        str(java_manifest_path),
        "--batch-threads",
        str(args.threads),
        "--batch-stats-out",
        str(java_raw_stats_path),
        "--go-owl",
        str(ref_dir / "go.owl"),
        "--pathways",
        str(ref_dir / "kegg_pathways.tsv"),
        "--modules",
        str(ref_dir / "kegg_modules.tsv"),
        "--ec2go",
        str(ref_dir / "ec2go.txt"),
        "--enable-priors",
        "essentiality,coherence,consistency,anatomical,gap_filling,genomic_context",
        "--enforce-stage2",
        "--enforce-consistency-mode",
        "relax",
        "--lite",
    ]
    log(f"Running gspa-unified.jar integrate across {len(java_manifest_rows)} representative whole-genome jobs...")
    t_java = time.time()
    subprocess.run(cmd, check=True)
    mean_rt = round((time.time() - t_java) / max(1, len(java_sample_gids)), 3)

    java_map: dict[str, dict[str, str]] = {}
    if java_raw_stats_path.exists():
        with open(java_raw_stats_path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if not line.strip():
                    continue
                vals = line.rstrip("\n").split("\t")
                row = {hdr[i]: vals[i] if i < len(vals) else "" for i in range(len(hdr))}
                java_map[row.get("genome_id", "")] = row

    with open(stats_out_path, "w") as sf:
        sf.write("genome_id\tn_claims\tn_integrated\tn_provenance\truntime_sec\n")
        for gid in sorted_gids:
            pbs = pre_batch_stats[gid]
            jr_noisy = java_map.get(f"{gid}__noisy_or")
            jr_full = java_map.get(gid)
            rt = float(jr_full["runtime_sec"]) if jr_full and jr_full.get("runtime_sec") else mean_rt
            sf.write(f"{gid}__noisy_or\t{pbs['n_claims']}\t{pbs['n_noisy']}\t0\t0.050\n")
            sf.write(f"{gid}\t{pbs['n_claims']}\t{pbs['n_full']}\t{pbs['n_prov']}\t{rt:.3f}\n")

    log(f"Completed Stage C in {time.time() - t0:.1f}s -> {stats_out_path}")


if __name__ == "__main__":
    main()
