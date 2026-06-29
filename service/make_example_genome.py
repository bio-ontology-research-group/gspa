#!/usr/bin/env python3
"""Build the example genomes shipped with DeepGOWeb's Genome (DeepGO-GSPA) tab:
per domain, a nucleotide FASTA + matching GFF3, assembled by reverse-translating
real proteins from a representative organism.

Reverse-translation (one fixed codon per residue, standard code) is deliberate:
gspa translates the CDS back to the exact same protein, so DeepGO-PlusPlus-Light's
DIAMOND step finds the real SwissProt/TrEMBL homolog and returns genuine GO
predictions -- while the build stays deterministic and vendors no third-party
sequence. Using standard codons (W->TGG) also sidesteps alternative genetic codes
(e.g. Mycoplasma's TGA=Trp): there is never an internal TGA to misread.

Four examples:
  * bacteria  -- Mycoplasma genitalium G37: the COMPLETE minimal bacterial
                 genome (~480 proteins)
  * phage     -- Enterobacteria phage T4: the COMPLETE viral genome (~290 proteins)
  * archaea   -- Methanocaldococcus jannaschii, ~150 best-annotated reviewed
                 proteins (the full ~1,800-protein proteome is too large for a
                 quick demo)
  * eukaryote -- Saccharomyces cerevisiae, ~150 best-annotated reviewed proteins
                 (full proteome ~6,000)

Bacteria + phage are complete genomes; archaea + eukaryote are representative
subsets. Accessions/sequences are fetched live from the UniProt REST API. Run
once (needs network); the emitted files are committed under
deepgoweb/apps/deepgo/examples/ and served by the example buttons.

Usage:
  make_example_genome.py                 # build all four
  make_example_genome.py bacteria phage
"""
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

CODON = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC', 'Q': 'CAG',
    'E': 'GAA', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA',
    'M': 'ATG', 'F': 'TTT', 'P': 'CCG', 'S': 'AGC', 'T': 'ACC', 'W': 'TGG',
    'Y': 'TAT', 'V': 'GTG',
}
STOP = 'TAA'
SPACER = 'TAGCTAGCATCGATCGTAGC'  # 20 bp intergenic filler (content irrelevant)
MAX_AA = 6000                    # keep essentially everything; drop only giants

# name -> dataset spec. complete=True streams the whole proteome; otherwise the
# `cap` best-annotated reviewed entries are taken.
# Bacteria uses E. coli K-12, a large representative subset rather than a complete
# genome: a complete E. coli is ~4,300 proteins (too slow for a quick demo), and
# the only feasible *complete* bacterium, Mycoplasma genitalium (~480), has such a
# reduced/divergent proteome that function-based inference mis-places it. A ~450
# E. coli subset is realistic AND infers Bacteria cleanly. Phage T4 stays complete.
DATASETS = {
    'bacteria':  dict(oid=83333,  contig='chromosome', complete=False, cap=450),  # E. coli K-12
    'phage':     dict(oid=10665,  contig='genome',     complete=True),   # phage T4 (complete)
    'archaea':   dict(oid=243232, contig='chromosome', complete=False, cap=150),  # M. jannaschii
    'eukaryote': dict(oid=559292, contig='chromosome', complete=False, cap=150),  # S. cerevisiae
}

HEADER_RE = re.compile(r'^>(?:\w+\|)?([A-Z0-9]+)(?:\|\S+)?\s')
GENE_RE = re.compile(r'\bGN=(\S+)')


def fetch_fasta(url):
    """Yield (accession, gene_or_accession, sequence) from a UniProt FASTA URL."""
    req = urllib.request.Request(url, headers={'Accept': 'text/plain'})
    with urllib.request.urlopen(req, timeout=120) as r:
        text = r.read().decode()
    acc, gene, seq = None, None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if acc and seq:
                yield acc, gene or acc, ''.join(seq)
            m = HEADER_RE.match(line + ' ')
            acc = m.group(1) if m else line[1:].split()[0]
            g = GENE_RE.search(line)
            gene = g.group(1) if g else None
            seq = []
        else:
            seq.append(line.strip())
    if acc and seq:
        yield acc, gene or acc, ''.join(seq)


def revtrans(seq):
    return ''.join(CODON[a] for a in seq) + STOP   # already starts with ATG (Met)


def proteome(spec):
    oid = spec['oid']
    if spec.get('complete'):
        q = urllib.parse.urlencode({'query': f'organism_id:{oid}', 'format': 'fasta'})
        url = f'https://rest.uniprot.org/uniprotkb/stream?{url_q(q)}'
    else:
        q = urllib.parse.urlencode({
            'query': f'organism_id:{oid} AND reviewed:true', 'format': 'fasta',
            'size': spec['cap'],
        })
        url = f'https://rest.uniprot.org/uniprotkb/search?{url_q(q)}'
    return list(fetch_fasta(url))


def url_q(q):
    # urlencode quotes spaces as '+'; the UniProt query parser is happy with that.
    return q


def build(name, spec):
    contig = spec['contig']
    raw = proteome(spec)
    seen, genes = set(), []
    skipped = 0
    for acc, gene, seq in raw:
        if not seq or any(a not in CODON for a in seq):
            skipped += 1
            continue
        if len(seq) > MAX_AA:
            skipped += 1
            continue
        base = gene if gene and gene not in seen else acc
        seen.add(base)
        genes.append((base, acc, revtrans(seq)))
    print(f'[{name}] {len(genes)} CDS kept, {skipped} skipped (non-standard residue / >{MAX_AA} aa)',
          file=sys.stderr)

    fasta, gff = [], ['##gff-version 3']
    seqparts, n = [], 0
    for gene, acc, nt in genes:
        seqparts.append(SPACER)
        start = sum(len(s) for s in seqparts) + 1          # 1-based
        seqparts.append(nt)
        end = start + len(nt) - 1
        n += 1
        attrs = f'ID={gene};locus_tag={contig}_{gene};product={gene}'
        gff.append(f'{contig}\tgspa-example\tgene\t{start}\t{end}\t.\t+\t.\tID=gene_{contig}_{n}')
        gff.append(f'{contig}\tgspa-example\tCDS\t{start}\t{end}\t.\t+\t0\t{attrs};Parent=gene_{contig}_{n}')
    seqparts.append(SPACER)
    full = ''.join(seqparts)
    fasta.append(f'>{contig}')
    fasta += [full[i:i + 70] for i in range(0, len(full), 70)]

    for prefix in ([name, 'example'] if name == 'bacteria' else [name]):
        with open(f'{prefix}_genome.fna', 'w') as fh:
            fh.write('\n'.join(fasta) + '\n')
        with open(f'{prefix}_genome.gff3', 'w') as fh:
            fh.write('\n'.join(gff) + '\n')
    print(f'  wrote {name}_genome.fna + .gff3 ({n} genes, {len(full)} bp)', file=sys.stderr)


def main():
    for name in (sys.argv[1:] or list(DATASETS)):
        if name not in DATASETS:
            print(f'unknown dataset: {name}', file=sys.stderr)
            continue
        print(f'== building {name} ==', file=sys.stderr)
        build(name, DATASETS[name])


if __name__ == '__main__':
    main()
