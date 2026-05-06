#!/usr/bin/env python3
"""Build a Swiss-Prot exclusion list covering every species in the panel.

For each species in panel_manifest.tsv, match on species name or taxon ID
in Swiss-Prot headers. Union of all matches is the exclude list we pass to
filter_fasta_by_exclude.py when building the leave-panel-out DIAMOND DB.

Header format we match:
  >sp|ACCESSION|NAME_SPECIES ... OS=Species name ... OX=taxid ...
"""
import argparse
import re
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--sprot-fasta', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    tax_ids = set()
    species_patterns = []
    with open(args.manifest) as f:
        next(f)
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            species = parts[1].strip()
            tax = parts[3].strip()
            tax_ids.add(tax)
            # Match the GENUS + species binomial from the common name, relaxed
            m = re.match(r'(\w+)\s+(\w+)', species)
            if m:
                species_patterns.append(f'{m.group(1)} {m.group(2)}')
    print(f'[info] {len(tax_ids)} taxon IDs, {len(species_patterns)} species patterns',
          file=sys.stderr)

    # Match either OS=<species_prefix> or OX=<taxid>
    ox_pat = re.compile(r'OX=(\d+)')
    os_pat_re = re.compile(r'OS=([^=]*?)(?:\s+[A-Z]{2}=|$)')

    excl = set()
    n_total = 0
    with open(args.sprot_fasta) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            n_total += 1
            header = line.rstrip('\n')
            m_ox = ox_pat.search(header)
            if m_ox and m_ox.group(1) in tax_ids:
                m_acc = re.match(r'>\w+\|([^|]+)\|', header)
                if m_acc:
                    excl.add(m_acc.group(1))
                continue
            m_os = os_pat_re.search(header)
            if m_os:
                os_name = m_os.group(1).strip()
                for sp in species_patterns:
                    if os_name.startswith(sp):
                        m_acc = re.match(r'>\w+\|([^|]+)\|', header)
                        if m_acc:
                            excl.add(m_acc.group(1))
                        break
    print(f'[info] {len(excl)} accessions excluded out of {n_total} Swiss-Prot entries',
          file=sys.stderr)

    with open(args.out, 'w') as f:
        for acc in sorted(excl):
            f.write(acc + '\n')
    print(f'[info] wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
