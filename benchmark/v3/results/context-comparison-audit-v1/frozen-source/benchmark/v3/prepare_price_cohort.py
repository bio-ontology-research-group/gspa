"""Prepare Price cohort family membership and metadata from sequence-similarity clusters.

Consumes the label-free layout directory (prepare_price_layout.py), an MMseqs2
easy-cluster membership TSV over the namespaced combined FASTA, and a frozen NCBI
taxonomy EFETCH XML. Produces the two files the production cohort engine reads:

  cohort_orthogroups.tsv  genome_id, protein_id, orthogroup_id, subfamily_id, weight
                          (CohortOrthogroupReader.readCohortFile)
  cohort_metadata.json    clade / strainCluster / assemblyStatus / contigs /
                          familySearchComplete / familySearchCoverage
                          (CohortMetadata.read, strict schema)

What the outputs are, and are not:
  * Families are sequence-similarity clusters under the recorded MMseqs2 protocol.
    They are not established orthologs; paralogs may share a family.
  * clade is the NCBI genus taxid and strainCluster the NCBI species taxid, both read
    from the frozen lineage. The species block is a conservative proxy for strain
    redundancy, not an ANI cluster. A rank that the lineage lacks stays absent
    (unknown); nothing is guessed.
  * familySearchComplete is true only when every acquired protein of the genome occurs
    exactly once in the cluster membership and no member ID is foreign. Absence of a
    family is detection-limited to this cluster protocol; it says nothing about remote
    homologs.
  * Only identifiers, assembly and contig facts enter the outputs. Functional GFF
    annotations and labels never do.

Protein IDs are exact-version IDs. The combined FASTA names a protein
``<assembly>__<protein_id>``; the namespace is split on the first ``__`` only. The same
protein ID in two genomes is two distinct proteins. Cluster IDs are stable hashes of
the representative's namespaced ID (no colon: the production reader strips ``ns:``).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA = 'gspa.price_cohort/1'
NAMESPACE_SEP = '__'
HEADER = ['genome_id', 'protein_id', 'orthogroup_id', 'subfamily_id', 'weight']

INTERPRETATION = (
    'Families are sequence-similarity clusters under the recorded MMseqs2 protocol, not established orthologs; '
    'paralogs can share a family. clade is the NCBI genus taxid and strainCluster the NCBI species taxid: the '
    'species block is a conservative redundancy proxy, not an ANI cluster, and strains/species are not independent '
    'units. familySearchComplete certifies only that every acquired protein appears exactly once in the cluster '
    'membership; absence of a family is detection-limited to this protocol and is not evidence of absence of remote '
    'homologs. No functional annotation or label was read.')


class CohortError(ValueError):
    pass


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def split_namespace(namespaced):
    """'<assembly>__<protein>' -> (assembly, protein), splitting on the first separator only."""
    genome, sep, protein = namespaced.partition(NAMESPACE_SEP)
    if not sep or not genome or not protein:
        raise CohortError(f'not a namespaced protein ID (<assembly>{NAMESPACE_SEP}<protein>): {namespaced!r}')
    return genome, protein


def stable_id(prefix, *parts):
    return f'{prefix}_{hashlib.sha256(chr(31).join(parts).encode()).hexdigest()[:12]}'


# ----------------------------------------------------------------------------- layout

def read_layout(layout_dir):
    """Verify the layout directory against its report and return acquired proteins per genome."""
    layout_dir = Path(layout_dir)
    report = json.loads((layout_dir / 'report.json').read_text())
    if report.get('schema') != 'gspa.price_layout/1':
        raise CohortError(f'unsupported layout report schema: {report.get("schema")!r}')
    metadata = json.loads((layout_dir / 'cohort_metadata.json').read_text())
    assemblies = report['assemblies']
    if set(metadata.get('genomes', {})) != set(assemblies):
        raise CohortError('layout cohort_metadata.json and report.json describe different assemblies')
    fasta = layout_dir / 'all_proteins.faa'
    if digest(fasta) != report['all_proteins_sha256']:
        raise CohortError('all_proteins.faa does not match the layout report digest')
    proteins, hashes = {}, {}
    for asm, rec in assemblies.items():
        seqfile = layout_dir / f'{asm}.sequence_hashes.tsv'
        if digest(seqfile) != rec['sequence_hashes_sha256']:
            raise CohortError(f'{seqfile.name} does not match the layout report digest')
        layout = layout_dir / f'{asm}.layout.tsv'
        if digest(layout) != rec['layout_sha256']:
            raise CohortError(f'{layout.name} does not match the layout report digest')
        rows = [line.rstrip('\n').split('\t') for line in seqfile.read_text().splitlines()]
        if rows[0] != ['protein_id', 'length', 'sequence_sha256']:
            raise CohortError(f'{seqfile.name}: unexpected header')
        ids = [r[0] for r in rows[1:]]
        if len(set(ids)) != len(ids):
            raise CohortError(f'{seqfile.name}: duplicate protein IDs')
        if len(ids) != rec['proteins']:
            raise CohortError(f'{seqfile.name}: {len(ids)} proteins, report says {rec["proteins"]}')
        proteins[asm] = set(ids)
        for r in rows[1:]:
            hashes[(asm, r[0])] = r[2]
        for line in layout.read_text().splitlines()[1:]:
            if line.split('\t')[0] not in proteins[asm]:
                raise CohortError(f'{layout.name}: coordinate row for a protein that was not acquired')
    audit_combined_fasta(fasta, proteins, hashes)
    return report, metadata, proteins


def audit_combined_fasta(path, proteins, hashes):
    """The combined FASTA must contain exactly the acquired namespaced proteins, with matching sequence digests."""
    seen, pid, parts = set(), None, []

    def finish():
        if pid is None:
            return
        asm, protein = split_namespace(pid)
        if asm not in proteins or protein not in proteins[asm]:
            raise CohortError(f'combined FASTA protein {pid!r} was not acquired')
        if pid in seen:
            raise CohortError(f'combined FASTA repeats {pid!r}')
        seen.add(pid)
        if hashlib.sha256(''.join(parts).encode()).hexdigest() != hashes[(asm, protein)]:
            raise CohortError(f'combined FASTA sequence of {pid!r} does not match its recorded digest')

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                finish()
                pid, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    finish()
    expected = sum(len(v) for v in proteins.values())
    if len(seen) != expected:
        raise CohortError(f'combined FASTA has {len(seen)} proteins, layout acquired {expected}')


# ----------------------------------------------------------------------------- clusters

def read_clusters(path):
    """MMseqs2 easy-cluster TSV (representative<TAB>member). Returns [(rep, member), ...] as written."""
    pairs = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            f = line.rstrip('\n').split('\t')
            if len(f) != 2 or not f[0] or not f[1]:
                raise CohortError(f'{path}:{n}: expected representative<TAB>member')
            pairs.append((f[0], f[1]))
    if not pairs:
        raise CohortError(f'{path} is empty')
    return pairs


def audit_membership(pairs, proteins):
    """Exact-once coverage audit. Returns (member -> representative, audit dict)."""
    expected = {f'{asm}{NAMESPACE_SEP}{p}' for asm, ps in proteins.items() for p in ps}
    owner, repeated, foreign = {}, [], []
    for rep, member in pairs:
        for pid in (rep, member):
            try:
                asm, protein = split_namespace(pid)
            except CohortError:
                foreign.append(pid)
                continue
            if asm not in proteins or protein not in proteins[asm]:
                foreign.append(pid)
        if member in owner:
            repeated.append(member)
        else:
            owner[member] = rep
    foreign = sorted(set(foreign))
    missing = sorted(expected - set(owner))
    for rep in sorted(set(owner.values())):
        if owner.get(rep) != rep:
            raise CohortError(f'representative {rep!r} is not listed as a member of its own cluster')
    audit = {'expected_proteins': len(expected), 'membership_rows': len(pairs), 'distinct_members': len(owner),
             'missing_proteins': len(missing), 'repeated_members': len(set(repeated)), 'foreign_ids': len(set(foreign)),
             'examples': {'missing': missing[:10], 'repeated': sorted(set(repeated))[:10], 'foreign': sorted(set(foreign))[:10]}}
    audit['exactly_once'] = not missing and not repeated and not foreign and len(owner) == len(pairs) == len(expected)
    audit['_missing'] = missing
    return owner, audit


def family_ids(owner, prefix):
    return {member: stable_id(prefix, rep) for member, rep in owner.items()}


def subfamily_assignment(owner30, owner70, policy):
    """Subfamily id per protein, validating that every 70% subcluster nests in one 30% family."""
    sub_families = defaultdict(set)
    for member, rep70 in owner70.items():
        sub_families[rep70].add(owner30[member])
    broken = sorted(rep for rep, fams in sub_families.items() if len(fams) > 1)
    warnings = []
    if broken and policy == 'fail':
        raise CohortError(f'{len(broken)} subclusters straddle several families (nesting violated), e.g. {broken[:3]}; '
                          'rerun with --subfamily-policy combined to use the explicit family+subcluster signature')
    if broken:
        warnings.append(f'{len(broken)} subclusters straddle several 30% families; subfamily_id is the combined '
                        'family+subcluster signature for ALL proteins, so subfamilies are strictly nested by construction')
    if policy == 'combined':
        return {m: stable_id('sub', owner30[m], r) for m, r in owner70.items()}, warnings, len(broken)
    return {m: stable_id('sub', r) for m, r in owner70.items()}, warnings, 0


# ----------------------------------------------------------------------------- taxonomy

def read_taxonomy(path, requested):
    """Parse EFETCH XML; top-level TaxId must be exactly a requested id. Returns {taxid: info}."""
    root = ET.parse(path).getroot()
    requested = {int(t) for t in requested}
    found = {}
    for taxon in root.findall('Taxon'):
        taxid = int(taxon.findtext('TaxId'))
        aka = {int(a.text) for a in taxon.findall('AkaTaxIds/TaxId')}
        if taxid not in requested:
            raise CohortError(f'taxonomy top-level TaxId {taxid} was not requested' +
                              (f' (merged from {sorted(aka & requested)})' if aka & requested else ''))
        if taxid in found:
            raise CohortError(f'taxonomy repeats TaxId {taxid}')
        ranks = {}
        for t in taxon.findall('LineageEx/Taxon'):
            rank = t.findtext('Rank')
            if rank in ('genus', 'species'):
                ranks.setdefault(rank, []).append(int(t.findtext('TaxId')))
        own = taxon.findtext('Rank')
        if own in ('genus', 'species') and taxid not in ranks.get(own, []):
            ranks.setdefault(own, []).append(taxid)
        # Lineage order is root-to-leaf; a rank repeated in one lineage is ambiguous and stays unknown.
        found[taxid] = {'name': taxon.findtext('ScientificName'), 'rank': own,
                        'genus': ranks['genus'][0] if len(ranks.get('genus', [])) == 1 else None,
                        'species': ranks['species'][0] if len(ranks.get('species', [])) == 1 else None}
    return found, sorted(requested - set(found))


def parse_mmseqs(log_path, version, command):
    if log_path:
        lines = Path(log_path).read_text().splitlines()
        version = next((l.split(':', 1)[1].strip() for l in lines if l.startswith('MMseqs Version')), None)
        command = next((l.strip() for l in lines if l.startswith(('easy-cluster ', 'cluster ', 'linclust '))), None)
        if command:
            command = 'mmseqs ' + command
    if not version or not command:
        raise CohortError('MMseqs2 version and command are required (--mmseqs-log, or --mmseqs-version and --mmseqs-command)')
    tokens = command.split()
    params = {t.lstrip('-'): v for t, v in zip(tokens, tokens[1:]) if t.startswith('-') and not v.startswith('-')}
    return {'version': version.strip(), 'command': command, 'parameters': params}


# ----------------------------------------------------------------------------- preparation

def prepare(layout_dir, clusters, taxonomy, out, mmseqs_log=None, mmseqs_version=None, mmseqs_command=None,
            subclusters=None, subfamily_policy='fail', allow_partial=False):
    out = Path(out)
    report_in, metadata_in, proteins = read_layout(layout_dir)
    protocol = parse_mmseqs(mmseqs_log, mmseqs_version, mmseqs_command)

    pairs = read_clusters(clusters)
    owner30, audit = audit_membership(pairs, proteins)
    if audit['foreign_ids'] or audit['repeated_members']:
        raise CohortError(f'cluster membership is invalid and cannot be certified: {audit["foreign_ids"]} foreign IDs, '
                          f'{audit["repeated_members"]} repeated members, e.g. {audit["examples"]}')
    if audit['missing_proteins'] and not allow_partial:
        raise CohortError(f'cluster membership covers {audit["distinct_members"]} of {audit["expected_proteins"]} acquired proteins '
                          f'({audit["missing_proteins"]} missing); refusing to certify (use --allow-partial for an uncertified output)')
    missing = set(audit.pop('_missing'))
    family = family_ids(owner30, 'fam30')

    sub, sub_warnings, sub_broken, sub_audit = {}, [], 0, None
    if subclusters:
        owner70, sub_audit = audit_membership(read_clusters(subclusters), proteins)
        sub_audit.pop('_missing')
        if sub_audit['foreign_ids'] or sub_audit['repeated_members'] or sub_audit['missing_proteins']:
            raise CohortError(f'subcluster membership is not exact-once over the acquired proteins: {sub_audit}')
        sub, sub_warnings, sub_broken = subfamily_assignment(owner30, owner70, subfamily_policy)

    taxids = {asm: rec.get('taxid') for asm, rec in report_in['assemblies'].items()}
    taxa, absent = read_taxonomy(taxonomy, {t for t in taxids.values() if t is not None})

    genomes, unresolved, per_genome = {}, {}, {}
    for asm in sorted(proteins):
        entry = dict(metadata_in['genomes'][asm])
        taxid = taxids[asm]
        info = taxa.get(taxid) if taxid is not None else None
        reasons = []
        if info is None:
            reasons.append('no taxid in layout report' if taxid is None else f'taxid {taxid} absent from frozen taxonomy')
        else:
            if info['genus'] is not None:
                entry['clade'] = f'ncbi_genus_{info["genus"]}'
            else:
                reasons.append('no unique genus rank in lineage')
            if info['species'] is not None:
                entry['strainCluster'] = f'ncbi_species_{info["species"]}'
            else:
                reasons.append('no unique species rank in lineage')
        if reasons:
            unresolved[asm] = reasons
        covered = len(proteins[asm]) - sum(1 for p in proteins[asm] if f'{asm}{NAMESPACE_SEP}{p}' in missing)
        coverage = covered / len(proteins[asm])
        entry['familySearchCoverage'] = round(coverage, 6) if coverage < 1 else 1.0
        entry['familySearchComplete'] = bool(audit['exactly_once']) and covered == len(proteins[asm])
        genomes[asm] = entry
        per_genome[asm] = {'proteins': len(proteins[asm]), 'in_clusters': covered, 'taxid': taxid}

    out.mkdir(parents=True, exist_ok=False)
    rows = 0
    with (out / 'cohort_orthogroups.tsv').open('w') as fh:
        fh.write('\t'.join(HEADER) + '\n')
        for asm in sorted(proteins):
            for protein in sorted(proteins[asm]):
                pid = f'{asm}{NAMESPACE_SEP}{protein}'
                if pid not in family:
                    continue
                fh.write(f'{asm}\t{protein}\t{family[pid]}\t{sub.get(pid, "-")}\t1\n')
                rows += 1
    (out / 'cohort_metadata.json').write_text(json.dumps({'genomes': genomes}, indent=2) + '\n')

    sizes = Counter(owner30.values())
    genomes_per_family, paralogous = defaultdict(set), set()
    per_genome_family = Counter()
    for member, rep in owner30.items():
        asm = split_namespace(member)[0]
        genomes_per_family[rep].add(asm)
        per_genome_family[(rep, asm)] += 1
    paralogous = {rep for (rep, _), n in per_genome_family.items() if n > 1}
    stats = {'families': len(sizes), 'singleton_families': sum(1 for n in sizes.values() if n == 1),
             'largest_family': max(sizes.values()),
             'families_in_more_than_one_genome': sum(1 for g in genomes_per_family.values() if len(g) > 1),
             'families_with_in_genome_paralogs': len(paralogous)}
    with (out / 'families.tsv').open('w') as fh:
        fh.write('orthogroup_id\trepresentative\tproteins\tgenomes\n')
        for rep, n in sorted(sizes.items()):
            fh.write(f'{stable_id("fam30", rep)}\t{rep}\t{n}\t{len(genomes_per_family[rep])}\n')

    result = {
        'schema': SCHEMA,
        'certified_family_search': bool(audit['exactly_once']),
        'family_unit': 'sequence-similarity cluster (not an established ortholog group)',
        'mmseqs': protocol,
        'cluster_audit': audit,
        'family_stats': stats,
        'subfamilies': ({'available': True, 'policy': subfamily_policy, 'nesting_violations': sub_broken, 'audit': sub_audit,
                         'warnings': sub_warnings, 'sha256': digest(subclusters)} if subclusters else
                        {'available': False, 'reason': 'no subcluster membership supplied; subfamily_id is "-"'}),
        'taxonomy': {'requested_taxids': sorted({t for t in taxids.values() if t is not None}), 'resolved': len(taxa),
                     'absent_from_xml': absent, 'unresolved_genomes': unresolved,
                     'clade': 'NCBI genus taxid from the frozen lineage', 'strainCluster': 'NCBI species taxid (proxy, not ANI)'},
        'genomes': per_genome,
        'output_rows': rows,
        'inputs': {'layout_report_sha256': digest(Path(layout_dir) / 'report.json'),
                   'layout_all_proteins_sha256': report_in['all_proteins_sha256'],
                   'layout_metadata_sha256': digest(Path(layout_dir) / 'cohort_metadata.json'),
                   'clusters_sha256': digest(clusters), 'taxonomy_sha256': digest(taxonomy),
                   'mmseqs_log_sha256': digest(mmseqs_log) if mmseqs_log else None},
        'outputs': {n: digest(out / n) for n in ('cohort_orthogroups.tsv', 'cohort_metadata.json', 'families.tsv')},
        'interpretation': INTERPRETATION,
    }
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--layout', required=True, help='directory written by prepare_price_layout.py')
    p.add_argument('--clusters', required=True, help='MMseqs2 easy-cluster <prefix>_cluster.tsv')
    p.add_argument('--taxonomy', required=True, help='frozen NCBI taxonomy EFETCH XML')
    p.add_argument('--out', required=True, help='new output directory')
    p.add_argument('--mmseqs-log', help='easy-cluster log (supplies version and command)')
    p.add_argument('--mmseqs-version')
    p.add_argument('--mmseqs-command')
    p.add_argument('--subclusters', help='optional 70%% subcluster membership TSV')
    p.add_argument('--subfamily-policy', choices=['fail', 'combined'], default='fail')
    p.add_argument('--allow-partial', action='store_true', help='write an uncertified output when proteins are missing')
    a = p.parse_args()
    r = prepare(a.layout, a.clusters, a.taxonomy, a.out, a.mmseqs_log, a.mmseqs_version, a.mmseqs_command,
                a.subclusters, a.subfamily_policy, a.allow_partial)
    print(json.dumps({'certified': r['certified_family_search'], 'rows': r['output_rows'], **r['family_stats']}))
