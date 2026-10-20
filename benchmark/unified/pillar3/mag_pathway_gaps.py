#!/usr/bin/env python3
"""
Pathway gap resolution across uncultured metagenome-assembled genomes (P14).

Regenerates the Suppl. Table S11 analysis ("pathway gap resolution across
uncultured MAGs") from first principles on the Empty Quarter high-quality
MAGs.  Standard library only (runs with any python3 >= 3.8 on ibex).

Subcommands
-----------
select     choose the MAG set (CheckM2 >= 90 / <= 5, completed gapseq run,
           GTDB classification, Bakta proteome + GFF3) and draw a seeded,
           phylum-stratified sample; writes mag_list.tsv + eligible table.
prep       per MAG: Bakta GFF3 -> genome layout + operons (locus_tag ids);
           gapseq draft vs final SBML -> gap reactions (reactions gapseq had
           to add to reach growth), gaps JSONL for `gspa integrate --gaps`,
           gapseq gene-level assignments mapped onto Bakta loci (control);
           MetaCyc pathway -> EC table for the gapseq-predicted pathways.
score      per MAG: count gap fills from Stage 0 / Stage 1+2A(+RLGC)
           posteriors, MetaCyc pathway completeness and essential-function
           completeness before/after, and the gapseq-gene control.
aggregate  per-MAG rows -> per-phylum table.

Definitions (also in paper-reproducibility/pipelines/P14_mag_pathway_gaps/README.md)
-----------
gap reaction   SEED reaction present in the gap-filled gapseq model
               (<mag>.xml) but absent from the draft (<mag>-draft.xml);
               exchange / demand / sink / biomass reactions excluded.
fillable gap   gap reaction with >= 1 complete EC number (no '-') that maps
               to >= 1 GO term through ec2go.
gap function   for a fillable gap: its complete EC numbers, their ec2go GO
               terms and all GO descendants of those terms (generic catalytic
               roots are never expanded nor matched) -- the same matching set
               GapFillingPrior uses.
fill @ q       (gap reaction, Bakta locus) pair such that the locus carries a
               gap function with posterior >= q in the Stage 1+2A output, or
               RLGC suggested the locus for the gap with q_value >= q.
new fill @ q   fill @ q whose locus had Stage 0 (noisy-OR only) posterior < 0.30
               for every gap function (not even a candidate before), so the
               q = 0.50 fills are a subset of the q = 0.30 fills.
               candidate_fills_q030 / hq_fills_q050 count new fills.
promoted       (gap, locus) pairs with Stage 0 posterior in [0.30, 0.50) that
               reach >= 0.50 after Stage 1+2A (reported separately).
gap resolved   fillable gap with >= 1 fill @ q and no Stage 0 locus >= q.
control        each new fill is classified by what gapseq itself had for the
               reaction: a gene in the gap-filled model's GPR that overlaps
               the GSPA locus (gapseq_gene_same_locus) or not
               (gapseq_gene_other_locus); no GPR (gapseq only added the
               reaction) but a sub-threshold gapseq tblastn hit
               (Reactions.tbl, status bad_blast etc.) overlapping the GSPA
               locus (gapseq_subthr_hit_same_locus) or elsewhere
               (gapseq_subthr_hit_other_locus); or no gapseq hit at all
               (gapseq_no_hit).
pathway        mean over gapseq-predicted MetaCyc pathways (Pathways.tbl
completeness   Prediction == true) of the fraction of the pathway's EC-bearing
               reactions (Reactions.tbl) that have a locus with posterior
               >= 0.5 for one of the reaction's functions (complete EC, its
               ec2go GO term or a descendant).  'before' = Stage 0; 'after' =
               Stage 1+2A plus RLGC suggestions with q >= 0.5.
essential      fraction of the bundled essential-function GO terms (bacteria
completeness   or archaea profile) carried, or via a descendant, by >= 1
               locus with posterior >= 0.5.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

GENERIC_ROOTS = {'GO:0003824', 'GO:0016491', 'GO:0016740', 'GO:0016787',
                 'GO:0016829', 'GO:0016853', 'GO:0016874', 'GO:0003674',
                 'GO:0008150', 'GO:0005575'}
EC_FULL = re.compile(r'^\d+\.\d+\.\d+\.\d+$')
QS = (0.30, 0.50)
CONTROL = ('gapseq_gene_same_locus', 'gapseq_gene_other_locus', 'gapseq_subthr_hit_same_locus', 'gapseq_subthr_hit_other_locus', 'gapseq_no_hit')


def xopen(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def read_tsv(path):
    with xopen(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, cols):
    with open(path, 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for r in rows:
            out.write('\t'.join(fmt(r.get(c, '')) for c in cols) + '\n')


def fmt(v):
    if isinstance(v, float):
        return 'NA' if math.isnan(v) else f'{v:.4f}'
    return str(v)


# ---------------------------------------------------------------------------
# select
# ---------------------------------------------------------------------------
def bakta_dir_for(fa_path, bin_id):
    # .../<sample>/wgs/binning/single_sample/output_bins/<bin>.fa.gz
    wgs = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(fa_path))))
    return os.path.join(wgs, 'annotation', 'single_sample', 'bakta', bin_id)


def cmd_select(a):
    ck = {}
    for r in read_tsv(a.checkm2):
        ck[f"{r['sample']}__{r['bin']}"] = r
    tax = {}
    for f in a.gtdb:
        for r in read_tsv(f):
            tax[r['user_genome']] = r['classification']
    rows, reasons = [], Counter()
    with open(a.hq_list) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            fa, sample, bin_id = f[0], f[1], f[2]
            mag = f'{sample}__{bin_id}'
            c = ck.get(mag)
            if c is None:
                reasons['no_checkm2'] += 1
                continue
            comp, cont = float(c['completeness']), float(c['contamination'])
            if comp < 90 or cont > 5:
                reasons['not_hq_checkm2'] += 1
                continue
            cl = tax.get(mag, '')
            ranks = dict(x.split('__', 1) for x in cl.split(';') if '__' in x)
            if not ranks.get('p'):
                reasons['no_gtdb_phylum'] += 1
                continue
            gd = os.path.join(a.gapseq_out, mag)
            need_g = [f'{mag}.xml', f'{mag}-draft.xml', f'{mag}-all-Pathways.tbl', f'{mag}-all-Reactions.tbl']
            if not all(os.path.isfile(os.path.join(gd, x)) and os.path.getsize(os.path.join(gd, x)) > 0 for x in need_g):
                reasons['gapseq_incomplete'] += 1
                continue
            bd = bakta_dir_for(fa, bin_id)
            if not all(os.path.isfile(os.path.join(bd, f'{bin_id}.{e}')) for e in ('faa', 'gff3', 'fna')):
                reasons['no_bakta'] += 1
                continue
            dom = ranks.get('d', '')
            rows.append({
                'mag_id': mag, 'sample': sample, 'bin': bin_id,
                'completeness': comp, 'contamination': cont,
                'domain': dom, 'phylum': ranks['p'], 'class': ranks.get('c', ''),
                'order': ranks.get('o', ''), 'family': ranks.get('f', ''),
                'genus': ranks.get('g', ''), 'species': ranks.get('s', ''),
                'gtdb_classification': cl,
                'taxon': 'NCBITaxon_2157' if dom == 'Archaea' else 'NCBITaxon_2',
                'essential_profile': 'archaea' if dom == 'Archaea' else 'bacteria',
                'mag_fasta': fa, 'bakta_dir': bd, 'gapseq_dir': gd,
            })
    cols = ['mag_id', 'sample', 'bin', 'completeness', 'contamination', 'domain', 'phylum',
            'class', 'order', 'family', 'genus', 'species', 'gtdb_classification', 'taxon',
            'essential_profile', 'mag_fasta', 'bakta_dir', 'gapseq_dir']
    write_tsv(a.eligible_out, rows, cols)
    by_ph = defaultdict(list)
    for r in rows:
        by_ph[r['phylum']].append(r)
    print(f'eligible MAGs: {len(rows)}; excluded: {dict(reasons)}')
    print('eligible per phylum:', {k: len(v) for k, v in sorted(by_ph.items(), key=lambda kv: -len(kv[1]))})
    if len(rows) <= a.full_set_max:
        chosen = rows
        print(f'<= {a.full_set_max} eligible: taking all')
    else:
        # Stratified allocation: phyla with >= min_phylum eligible MAGs each get
        # `base` slots; the remaining slots go proportionally to phylum size
        # (largest-remainder); MAGs drawn uniformly within phylum (seeded).
        strata = {p: v for p, v in by_ph.items() if len(v) >= a.min_phylum}
        alloc = {p: min(a.base, len(v)) for p, v in strata.items()}
        left = a.n - sum(alloc.values())
        if left < 0:
            raise SystemExit('base * n_strata exceeds n')
        tot = sum(len(v) - alloc[p] for p, v in strata.items())
        quotas = {p: left * (len(v) - alloc[p]) / tot for p, v in strata.items()}
        add = {p: int(math.floor(q)) for p, q in quotas.items()}
        rem = left - sum(add.values())
        for p in sorted(quotas, key=lambda p: (-(quotas[p] - add[p]), p))[:rem]:
            add[p] += 1
        rng = random.Random(a.seed)
        chosen = []
        for p in sorted(strata):
            k = min(len(strata[p]), alloc[p] + add[p])
            pool = sorted(strata[p], key=lambda r: r['mag_id'])
            chosen.extend(rng.sample(pool, k))
        print('allocation:', {p: alloc[p] + add[p] for p in sorted(strata)})
    chosen.sort(key=lambda r: (r['phylum'], r['mag_id']))
    write_tsv(a.out, chosen, cols)
    print(f'selected {len(chosen)} MAGs across {len({r["phylum"] for r in chosen})} phyla -> {a.out}')


# ---------------------------------------------------------------------------
# prep
# ---------------------------------------------------------------------------
def load_ec2go(path):
    out = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            m = re.match(r'^EC:(\S+)\s*>\s*GO:.*;\s*(GO:\d{7})', line)
            if m:
                out[m.group(1)].add(m.group(2))
    return out


def parse_sbml(path):
    """Return {seed_rxn: {'ecs': set, 'genes': [gene ids]}} for a gapseq SBML model.
    Namespace-agnostic (gapseq writes SBML L3V1 or L3V2 core)."""
    rdf = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
    out = {}
    root = ET.parse(path).getroot()
    for rx in root.iter():
        if not rx.tag.endswith('}reaction'):
            continue
        rid = rx.get('id', '')
        m = re.match(r'^R_(rxn\d+)', rid)
        if not m:
            continue
        seed = m.group(1)
        ecs, genes = set(), []
        for el in rx.iter():
            if el.tag.endswith('}li'):
                res = el.get('{%s}resource' % rdf, '')
                if 'ec-code:' in res:
                    ecs.add(res.split('ec-code:')[1])
            elif el.tag.endswith('}geneProductRef'):
                g = next((v for k, v in el.attrib.items() if k.endswith('geneProduct')), None)
                if g:
                    genes.append(g)
        e = out.setdefault(seed, {'ecs': set(), 'genes': [], 'sbml_id': rid})
        e['ecs'] |= ecs
        e['genes'].extend(genes)
    return out


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with xopen(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf).upper()
    return seqs


def norm(s):
    return re.sub(r'[^A-Za-z0-9]', '_', s)


def read_tbl(path):
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if header is None:
                header = f
                continue
            f += [''] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows


def split_ecs(s):
    return [e for e in re.split(r'[/,;\s]+', s or '') if e and e not in ('NA', '-')]


def cmd_prep(a):
    os.makedirs(a.out_dir, exist_ok=True)
    mag, bin_id = a.mag_id, a.bin
    P = lambda name: os.path.join(a.out_dir, f'{mag}_{name}')
    # --- layout + operons from Bakta GFF3 (protein id = locus_tag = .faa id)
    faa_ids = set()
    with open(os.path.join(a.bakta_dir, f'{bin_id}.faa')) as fh:
        for line in fh:
            if line.startswith('>'):
                faa_ids.add(line[1:].split()[0])
    cds = []
    with open(os.path.join(a.bakta_dir, f'{bin_id}.gff3')) as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            m = re.search(r'(?:^|;)ID=([^;]+)', f[8])
            if not m or m.group(1) not in faa_ids:
                continue
            cds.append((f[0], int(f[3]), int(f[4]), f[6], m.group(1)))
    cds.sort(key=lambda r: (r[0], r[1]))
    with open(P('layout.tsv'), 'w') as out:
        out.write('protein_id\tcontig\tstart\tend\tstrand\n')
        for c, s, e, st, pid in cds:
            out.write(f'{pid}\t{c}\t{s}\t{e}\t{st}\n')
    operons, cur, prev = [], [], None
    for c, s, e, st, pid in cds:
        if prev and prev[0] == c and prev[3] == st and s - prev[2] <= 300:
            cur.append(pid)
        else:
            if len(cur) >= 2:
                operons.append(cur)
            cur = [pid]
        prev = (c, s, e, st)
    if len(cur) >= 2:
        operons.append(cur)
    with open(P('operons.tsv'), 'w') as out:
        for op in operons:
            out.write('\t'.join(op) + '\n')

    # --- contig map: original MAG contig name -> Bakta contig (by sequence)
    orig = read_fasta(a.mag_fasta)
    bak = read_fasta(os.path.join(a.bakta_dir, f'{bin_id}.fna'))
    by_seq = {hashlib.md5(s.encode()).hexdigest(): n for n, s in bak.items()}
    cmap = {}
    for n, s in orig.items():
        b = by_seq.get(hashlib.md5(s.encode()).hexdigest())
        if b:
            cmap[norm(n)] = b
    loci_by_contig = defaultdict(list)
    for c, s, e, st, pid in cds:
        loci_by_contig[c].append((s, e, pid))

    def gene_to_loci(gid):
        # gapsmith gene products are protein ids (G_<protein_id>, DIAMOND on the proteome)
        if gid and gid.startswith('G_') and gid[2:] in faa_ids:
            return 'protein', [gid[2:]]
        # gapseq gene products: G_<contig>_<sstart>_<send> (tblastn coordinates)
        m = re.match(r'^G_(.+)_(\d+)_(\d+)$', gid or '')
        if not m:
            return None, []
        contig = cmap.get(norm(m.group(1)))
        if contig is None:
            return None, []
        lo, hi = sorted((int(m.group(2)), int(m.group(3))))
        hits = [pid for s, e, pid in loci_by_contig[contig] if s <= hi and e >= lo]
        return contig, hits

    # --- gap reactions = final model minus draft model
    final = parse_sbml(os.path.join(a.gapseq_dir, f'{mag}.xml'))
    draft = parse_sbml(os.path.join(a.gapseq_dir, f'{mag}-draft.xml'))
    ec2go = load_ec2go(a.ec2go)
    # MetaCyc pathway of a SEED reaction (via the dbhit column of Reactions.tbl)
    rtbl = read_tbl(os.path.join(a.gapseq_dir, f'{mag}-all-Reactions.tbl'))
    seed2pwy = defaultdict(set)
    seed2hits = defaultdict(list)   # gapseq tblastn hits (any status) of the reaction's reference proteins
    seed2ec_tbl = defaultdict(set)  # EC numbers of the MetaCyc reactions a SEED reaction is a hit of
    for r in rtbl:
        pw = (r.get('pathway') or '').strip('|')
        for s in (r.get('dbhit') or '').split():
            if not s.startswith('rxn'):
                continue
            if pw:
                seed2pwy[s].add(pw)
            seed2ec_tbl[s].update(split_ecs(r.get('ec')))
            if r.get('stitle') not in (None, '', 'NA') and (r.get('sstart') or '').isdigit() and (r.get('send') or '').isdigit():
                seed2hits[s].append(r)
    # Optional EC sources beyond the SBML annotations (gapsmith models carry none): the Reactions.tbl
    # column `ec` via `dbhit`, and the ModelSEED Enzyme Class alias table (rxn ids in columns 1-2, EC in 3).
    seed2ec_alias = defaultdict(set)
    if getattr(a, 'reaction_ec_aliases', None):
        with open(a.reaction_ec_aliases) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 3 or not EC_FULL.match(f[2]):
                    continue
                for s in re.split(r'[|;, ]+', f[0] + '|' + f[1]):
                    if s.startswith('rxn'):
                        seed2ec_alias[s].add(f[2])
    use_tbl_ecs = bool(getattr(a, 'ec_from_reactions_tbl', False))
    rank = {'good_blast': 3, 'bad_blast': 2}

    def hit_loci(seed):
        loci, best, bits = set(), '', 0.0
        for r in seed2hits.get(seed, ()):
            st = r.get('status', '')
            if rank.get(st, 1) > rank.get(best, 0):
                best = st
            try:
                bits = max(bits, float(r.get('bitscore') or 0))
            except ValueError:
                pass
            pid = r['stitle'].split()[0]
            if pid in faa_ids:          # gapsmith: subject is the protein itself
                loci.add(pid)
                continue
            contig = cmap.get(norm(r['stitle']))
            if contig is None:
                continue
            lo, hi = sorted((int(r['sstart']), int(r['send'])))
            loci.update(pid for s0, e0, pid in loci_by_contig[contig] if s0 <= hi and e0 >= lo)
        return loci, best, bits
    gaps_rows, n_jsonl = [], 0
    unmapped_genes = 0
    with open(P('gaps.jsonl'), 'w') as gj:
        for seed in sorted(set(final) - set(draft)):
            info = final[seed]
            ecs = set(info['ecs']) | seed2ec_alias.get(seed, set())
            if use_tbl_ecs:
                ecs |= seed2ec_tbl.get(seed, set())
            ecs = sorted(ecs)
            full = [e for e in ecs if EC_FULL.match(e)]
            gos = sorted({g for e in full for g in ec2go.get(e, ())})
            pw = sorted(seed2pwy.get(seed, ()))
            pathway_id = f'MetaCyc:{pw[0]}' if pw else f'biomass_demand:{seed}'
            genes = sorted(set(info['genes']))
            gl = set()
            for g in genes:
                c, hits = gene_to_loci(g)
                if c is None:
                    unmapped_genes += 1
                gl.update(hits)
            hl, hst, hbits = hit_loci(seed)
            gaps_rows.append({
                'reaction_id': seed, 'ecs': ','.join(ecs), 'complete_ecs': ','.join(full),
                'go_terms': ','.join(gos), 'pathways': ','.join(pw),
                'fillable': int(bool(full and gos)),
                'gapseq_has_gene': int(bool(genes)), 'gapseq_genes': ','.join(genes),
                'gapseq_loci': ','.join(sorted(gl)),
                'gapseq_hit_status': hst, 'gapseq_hit_bitscore': hbits,
                'gapseq_hit_loci': ','.join(sorted(hl)),
            })
            if not full:
                continue
            # one record per (reaction, EC, GO); GO may be null (prior matches the EC)
            for e in full:
                for g in (sorted(ec2go.get(e, ())) or [None]):
                    gj.write(json.dumps({'pathway_id': pathway_id, 'reaction_id': seed,
                                         'ec_number': e, 'go_term': g,
                                         'gapseq_guessed': bool(genes)}) + '\n')
                    n_jsonl += 1
    write_tsv(P('gap_reactions.tsv'), gaps_rows,
              ['reaction_id', 'ecs', 'complete_ecs', 'go_terms', 'pathways', 'fillable',
               'gapseq_has_gene', 'gapseq_genes', 'gapseq_loci', 'gapseq_hit_status',
               'gapseq_hit_bitscore', 'gapseq_hit_loci'])

    # --- predicted MetaCyc pathways and their EC-bearing reactions
    ptbl = read_tbl(os.path.join(a.gapseq_dir, f'{mag}-all-Pathways.tbl'))
    pred = {}
    for r in ptbl:
        if (r.get('Prediction') or '').strip().lower() == 'true':
            try:
                pred[r['ID'].strip('|')] = float(r.get('Completeness') or 'nan')
            except ValueError:
                pred[r['ID'].strip('|')] = float('nan')
    prx = defaultdict(lambda: defaultdict(set))
    for r in rtbl:
        pw = (r.get('pathway') or '').strip('|')
        if pw in pred:
            for e in split_ecs(r.get('ec')):
                prx[pw][r['rxn']].add(e)
    with open(P('pathway_reactions.tsv'), 'w') as out:
        out.write('pathway_id\tgapseq_completeness\treaction_id\tecs\n')
        for pw in sorted(pred):
            for rx in sorted(prx[pw]):
                out.write(f'{pw}\t{pred[pw]}\t{rx}\t{",".join(sorted(prx[pw][rx]))}\n')
    n_fill = sum(r['fillable'] for r in gaps_rows)
    n_gene = sum(r['gapseq_has_gene'] for r in gaps_rows)
    n_hit = sum(bool(r['gapseq_hit_loci']) for r in gaps_rows if r['fillable'])
    summ = {'mag_id': mag, 'n_proteins': len(faa_ids), 'n_cds_layout': len(cds),
            'n_operons': len(operons), 'n_contigs': len(bak), 'n_contigs_mapped': len(cmap),
            'n_final_rxns': len(final), 'n_draft_rxns': len(draft),
            'n_gap_rxns': len(gaps_rows), 'n_gap_fillable': n_fill,
            'n_gap_with_gapseq_gene': n_gene, 'n_fillable_gap_with_gapseq_hit_locus': n_hit, 'n_gapseq_genes_unmapped': unmapped_genes,
            'n_gap_jsonl': n_jsonl, 'n_pred_pathways': len(pred)}
    with open(P('prep_summary.json'), 'w') as out:
        json.dump(summ, out, indent=1)
    print(json.dumps(summ))


# ---------------------------------------------------------------------------
# score
# ---------------------------------------------------------------------------
def load_desc(go_dag):
    """go-dag.tsv rows: term <TAB> ancestor (reflexive closure)."""
    desc = defaultdict(set)
    with open(go_dag) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2:
                desc[f[1]].add(f[0])
    return desc


def load_post(path):
    """function_id -> {protein: posterior} (max over duplicate rows)."""
    by_term = defaultdict(dict)
    if not path or not os.path.exists(path):
        return by_term
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ip, it, ipp = header.index('protein_id'), header.index('function_id'), header.index('posterior_prob')
        for line in fh:
            f = line.rstrip('\n').split('\t')
            try:
                p = float(f[ipp])
            except (ValueError, IndexError):
                continue
            d = by_term[f[it]]
            if p > d.get(f[ip], -1):
                d[f[ip]] = p
    return by_term


def load_sugg(path):
    """reaction_id -> {protein: max q}."""
    out = defaultdict(dict)
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    for r in read_tsv(path):
        for pid, q in zip((r.get('protein_ids') or '').split(','), (r.get('q_values') or '').split(',')):
            try:
                qv = float(q)
            except ValueError:
                continue
            if pid and qv > out[r['reaction_id']].get(pid, -1):
                out[r['reaction_id']][pid] = qv
    return out


def func_set(ecs, ec2go, desc):
    terms = set()
    for e in ecs:
        if not EC_FULL.match(e):
            continue
        terms.add(e)
        for g in ec2go.get(e, ()):
            if g in GENERIC_ROOTS:
                continue
            terms.add(g)
            terms |= {d for d in desc.get(g, ()) if d not in GENERIC_ROOTS}
    return terms


def loci_at(terms, post, q):
    out = {}
    for t in terms:
        for pid, p in post.get(t, {}).items():
            if p >= q and p > out.get(pid, -1):
                out[pid] = p
    return out


def cmd_score(a):
    mag = a.mag_id
    P = lambda name: os.path.join(a.prep_dir, f'{mag}_{name}')
    ec2go = load_ec2go(a.ec2go)
    desc = load_desc(a.go_dag)
    s0 = load_post(a.stage0)
    s12 = load_post(a.stage12)
    sug = load_sugg(a.suggestions)
    gaps = read_tsv(P('gap_reactions.tsv'))
    fills_out = []
    res = {'mag_id': mag}
    for q in QS:
        tag = f'q{int(round(q * 100)):03d}'
        n_new = n_prior = n_rlgc = n_s0 = n_prom = 0
        resolved = covered0 = 0
        ctl = Counter()
        for g in gaps:
            if g['fillable'] != '1':
                continue
            terms = func_set(g['complete_ecs'].split(','), ec2go, desc)
            f0 = loci_at(terms, s0, q)
            f12 = loci_at(terms, s12, q)
            fr = {p: v for p, v in sug.get(g['reaction_id'], {}).items() if v >= q}
            after = set(f12) | set(fr)
            # baseline for a 'new' pair is fixed at Stage 0 >= 0.30 (not even a candidate
            # before), so fills at 0.50 are a subset of fills at 0.30
            f0_cand = loci_at(terms, s0, QS[0])
            new = sorted(after - set(f0_cand))
            n_s0 += len(f0)
            covered0 += int(bool(f0))
            if after and not f0:
                resolved += 1
            if q == 0.50:
                n_prom += sum(1 for pid in f12 if pid in f0_cand and pid not in f0)
            gl = set(filter(None, g['gapseq_loci'].split(',')))
            hl = set(filter(None, g['gapseq_hit_loci'].split(',')))
            for pid in new:
                n_new += 1
                src = '+'.join(x for x, ok in (('prior', pid in f12), ('rlgc', pid in fr)) if ok)
                n_prior += int(pid in f12)
                n_rlgc += int(pid in fr)
                if g['gapseq_has_gene'] == '1':
                    c = 'gapseq_gene_same_locus' if pid in gl else 'gapseq_gene_other_locus'
                elif hl:
                    c = 'gapseq_subthr_hit_same_locus' if pid in hl else 'gapseq_subthr_hit_other_locus'
                else:
                    c = 'gapseq_no_hit'
                ctl[c] += 1
                fills_out.append({'mag_id': mag, 'q': q, 'reaction_id': g['reaction_id'],
                                  'complete_ecs': g['complete_ecs'], 'pathways': g['pathways'],
                                  'locus': pid, 'posterior_stage12': f12.get(pid, ''),
                                  'rlgc_q': fr.get(pid, ''), 'source': src, 'control': c,
                                  'gapseq_loci': g['gapseq_loci'], 'gapseq_hit_loci': g['gapseq_hit_loci'],
                                  'gapseq_hit_status': g['gapseq_hit_status']})
        res[f'new_fills_{tag}'] = n_new
        res[f'new_fills_prior_{tag}'] = n_prior
        res[f'new_fills_rlgc_{tag}'] = n_rlgc
        res[f'stage0_fills_{tag}'] = n_s0
        res[f'gaps_covered_stage0_{tag}'] = covered0
        res[f'gaps_resolved_{tag}'] = resolved
        if q == 0.50:
            res['promoted_fills_q050'] = n_prom
        for c in CONTROL:
            res[f'fills_{c}_{tag}'] = ctl[c]
    res['n_gap_rxns'] = len(gaps)
    res['n_gap_fillable'] = sum(g['fillable'] == '1' for g in gaps)
    res['n_gap_fillable_gapseq_gene'] = sum(g['fillable'] == '1' and g['gapseq_has_gene'] == '1' for g in gaps)
    res['n_gap_fillable_gapseq_hit'] = sum(g['fillable'] == '1' and bool(g['gapseq_hit_loci']) for g in gaps)
    res['n_rlgc_suggested_rxns'] = len(sug)
    # RLGC diagnostic: does the top-ranked RLGC locus coincide with gapseq's sub-threshold hit?
    top1 = top1_hit = top_q = 0
    for g in gaps:
        cand = sug.get(g['reaction_id'])
        hl = set(filter(None, g['gapseq_hit_loci'].split(',')))
        if g['fillable'] != '1' or not cand:
            continue
        pid, qv = max(cand.items(), key=lambda kv: kv[1])
        top1 += 1
        top1_hit += int(pid in hl)
        top_q = max(top_q, qv)
    res['rlgc_top1_n'] = top1
    res['rlgc_top1_is_gapseq_hit_locus'] = top1_hit
    res['rlgc_max_q'] = float(top_q)

    # pathway completeness (q = 0.5)
    pw = defaultdict(dict)
    gs_comp = {}
    for r in read_tsv(P('pathway_reactions.tsv')):
        gs_comp[r['pathway_id']] = r['gapseq_completeness']
        ecs = [e for e in r['ecs'].split(',') if e]
        if any(EC_FULL.match(e) for e in ecs):
            pw[r['pathway_id']][r['reaction_id']] = ecs
    # RLGC suggestions >= 0.5 lift the gap reaction's ECs
    rlgc_ecs = set()
    for g in gaps:
        if any(v >= 0.5 for v in sug.get(g['reaction_id'], {}).values()):
            rlgc_ecs.update(e for e in g['complete_ecs'].split(',') if e)
    cache = {}

    def covered(ecs, post, extra=()):
        key = (tuple(ecs), id(post), bool(extra))
        if key not in cache:
            terms = func_set(ecs, ec2go, desc)
            cache[key] = bool(loci_at(terms, post, 0.5)) or any(e in extra for e in ecs)
        return cache[key]

    before, after = [], []
    for p, rx in pw.items():
        n = len(rx)
        before.append(sum(covered(e, s0) for e in rx.values()) / n)
        after.append(sum(covered(e, s12, rlgc_ecs) for e in rx.values()) / n)
    res['n_pred_pathways_scored'] = len(pw)
    res['pathway_completeness_before'] = sum(before) / len(before) if before else float('nan')
    res['pathway_completeness_after'] = sum(after) / len(after) if after else float('nan')
    gsv = [float(v) for v in gs_comp.values() if v not in ('', 'nan', 'NA')]
    res['gapseq_pathway_completeness'] = (sum(gsv) / len(gsv) / 100.0) if gsv else float('nan')

    # essential completeness (q = 0.5)
    ess = [l.split('\t')[0] for l in open(a.essential) if l.startswith('GO:')]

    def ess_cov(post):
        k = 0
        for t in ess:
            terms = ({t} | desc.get(t, set()))
            k += bool(loci_at(terms, post, 0.5))
        return k / len(ess) if ess else float('nan')
    res['n_essential_terms'] = len(ess)
    res['essential_completeness_before'] = ess_cov(s0)
    res['essential_completeness_after'] = ess_cov(s12)
    res['n_annot_stage0_q050'] = sum(1 for d in s0.values() for p in d.values() if p >= 0.5)
    res['n_annot_stage12_q050'] = sum(1 for d in s12.values() for p in d.values() if p >= 0.5)
    with open(a.out, 'w') as out:
        out.write('\t'.join(res) + '\n')
        out.write('\t'.join(fmt(v) for v in res.values()) + '\n')
    write_tsv(a.fills_out, fills_out, ['mag_id', 'q', 'reaction_id', 'complete_ecs', 'pathways', 'locus',
                                       'posterior_stage12', 'rlgc_q', 'source', 'control', 'gapseq_loci',
                                       'gapseq_hit_loci', 'gapseq_hit_status'])
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()}))


# ---------------------------------------------------------------------------
# aggregate
# ---------------------------------------------------------------------------
def cmd_aggregate(a):
    mags = {r['mag_id']: r for r in read_tsv(a.mag_list)}
    per = []
    missing = []
    for m, info in mags.items():
        f = os.path.join(a.work, m, f'{m}_score.tsv')
        if not os.path.exists(f):
            missing.append(m)
            continue
        r = read_tsv(f)[0]
        r.update({'domain': info['domain'], 'phylum': info['phylum'],
                  'completeness': info['completeness'], 'contamination': info['contamination'],
                  'gtdb_classification': info['gtdb_classification']})
        per.append(r)
    if missing:
        print(f'WARNING: {len(missing)} MAGs without score: {missing}', file=sys.stderr)
    ints = [k for k in per[0] if k.startswith(('new_fills', 'stage0_fills', 'gaps_', 'fills_', 'n_', 'promoted', 'rlgc_top1'))]
    floats = ['rlgc_max_q', 'pathway_completeness_before', 'pathway_completeness_after', 'gapseq_pathway_completeness',
              'essential_completeness_before', 'essential_completeness_after']
    lead = ['mag_id', 'domain', 'phylum', 'completeness', 'contamination', 'gtdb_classification']
    write_tsv(a.per_mag_out, per, lead + ints + floats)

    def mean(xs):
        xs = [x for x in xs if not math.isnan(x)]
        return sum(xs) / len(xs) if xs else float('nan')

    groups = defaultdict(list)
    for r in per:
        groups[(r['phylum'], r['domain'])].append(r)
    rows = []
    for key in sorted(groups, key=lambda k: (-len(groups[k]), k)) + [('ALL', 'ALL')]:
        rs = per if key[0] == 'ALL' else groups[key]
        S = lambda c: sum(int(float(r[c])) for r in rs)
        F = lambda c: mean([float(r[c]) for r in rs])
        n = len(rs)
        row = {'phylum_or_group': key[0], 'domain': key[1], 'n_mags': n,
               'candidate_fills_q030': S('new_fills_q030'), 'hq_fills_q050': S('new_fills_q050'),
               'mean_fills_per_mag_q030': S('new_fills_q030') / n,
               'mean_fills_per_mag_q050': S('new_fills_q050') / n,
               'pathway_completeness_before': F('pathway_completeness_before'),
               'pathway_completeness_after': F('pathway_completeness_after'),
               'essential_completeness_before': F('essential_completeness_before'),
               'essential_completeness_after': F('essential_completeness_after'),
               'gapseq_pathway_completeness': F('gapseq_pathway_completeness'),
               'mean_gap_rxns_per_mag': S('n_gap_rxns') / n,
               'mean_fillable_gaps_per_mag': S('n_gap_fillable') / n,
               'gaps_resolved_q030': S('gaps_resolved_q030'), 'gaps_resolved_q050': S('gaps_resolved_q050'),
               'gaps_covered_stage0_q050': S('gaps_covered_stage0_q050'),
               'fills_rlgc_q030': S('new_fills_rlgc_q030'), 'fills_rlgc_q050': S('new_fills_rlgc_q050'),
               'fills_prior_q030': S('new_fills_prior_q030'), 'fills_prior_q050': S('new_fills_prior_q050'),
               'promoted_fills_q050': S('promoted_fills_q050'),
               'rlgc_top1_n': S('rlgc_top1_n'), 'rlgc_top1_is_gapseq_hit_locus': S('rlgc_top1_is_gapseq_hit_locus'),
               'rlgc_max_q': max(float(r['rlgc_max_q']) for r in rs)}
        for tag in ('q030', 'q050'):
            for c in CONTROL:
                row[f'{c}_{tag}'] = S(f'fills_{c}_{tag}')
        rows.append(row)
    write_tsv(a.out, rows, list(rows[0].keys()))
    print(f'{len(per)} MAGs, {len(groups)} phyla -> {a.out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('select')
    s.add_argument('--hq-list', required=True)
    s.add_argument('--checkm2', required=True)
    s.add_argument('--gtdb', nargs='+', required=True)
    s.add_argument('--gapseq-out', required=True)
    s.add_argument('--n', type=int, default=100)
    s.add_argument('--full-set-max', type=int, default=150)
    s.add_argument('--min-phylum', type=int, default=10)
    s.add_argument('--base', type=int, default=4)
    s.add_argument('--seed', type=int, default=20261019)
    s.add_argument('--eligible-out', required=True)
    s.add_argument('--out', required=True)
    s.set_defaults(fn=cmd_select)
    p = sp.add_parser('prep')
    for k in ('--mag-id', '--bin', '--mag-fasta', '--bakta-dir', '--gapseq-dir', '--ec2go', '--out-dir'):
        p.add_argument(k, required=True)
    p.add_argument('--reaction-ec-aliases', default=None,
                   help='ModelSEED Enzyme Class alias table; adds its EC numbers to the gap reactions (for models without SBML EC annotations, e.g. gapsmith)')
    p.add_argument('--ec-from-reactions-tbl', action='store_true',
                   help='also take EC numbers from the Reactions.tbl rows whose dbhit names the gap reaction')
    p.set_defaults(fn=cmd_prep)
    c = sp.add_parser('score')
    for k in ('--mag-id', '--prep-dir', '--stage0', '--stage12', '--ec2go', '--go-dag', '--essential',
              '--out', '--fills-out'):
        c.add_argument(k, required=True)
    c.add_argument('--suggestions', default=None)
    c.set_defaults(fn=cmd_score)
    g = sp.add_parser('aggregate')
    for k in ('--mag-list', '--work', '--out', '--per-mag-out'):
        g.add_argument(k, required=True)
    g.set_defaults(fn=cmd_aggregate)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
