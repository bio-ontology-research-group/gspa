#!/usr/bin/env python3
"""Plan step S1: verify that the pipeline's ec2go file agrees with the EC cross-references of the pinned go.obo.

For every EC in ec2go: is it an `xref: EC:` of the same GO term in the obo? Lists ECs only in ec2go, only in the obo,
and ECs whose GO term differs. Optionally checks a list of ECs (one per line, e.g. the truth ECs) in both. Reads no
label; the EC list argument is a list of identifiers.
"""
import argparse, hashlib, json, re
from pathlib import Path


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_ec2go(path):
    m = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('EC:') and ';' in line:
            m.setdefault(line.split()[0][3:], set()).add(line.rsplit(';', 1)[1].strip())
    return m


def read_obo(path):
    m, cur, obsolete = {}, None, set()
    for line in open(path):
        line = line.rstrip('\n')
        if line == '[Term]':
            cur = None
        elif line.startswith('id: GO:'):
            cur = line[4:]
        elif line.startswith('xref: EC:') and cur:
            m.setdefault(line.split()[1][3:], set()).add(cur)
        elif line == 'is_obsolete: true' and cur:
            obsolete.add(cur)
    return m, obsolete


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--obo', required=True)
    ap.add_argument('--ec-list')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    e = read_ec2go(a.ec2go)
    o, obsolete = read_obo(a.obo)
    only_ec2go = sorted(set(e) - set(o))
    only_obo = sorted(set(o) - set(e))
    differ = sorted(ec for ec in set(e) & set(o) if e[ec] != o[ec])
    rep = {'schema': 'gspa.row11-ec2go-obo-check/1', 'ec2go': {'path': a.ec2go, 'sha256': sha(a.ec2go), 'ecs': len(e)},
           'obo': {'path': a.obo, 'sha256': sha(a.obo), 'ecs': len(o), 'terms_with_ec_xref': len({t for s in o.values() for t in s}),
                   'obsolete_terms_with_ec_xref': len({t for s in o.values() for t in s} & obsolete)},
           'only_in_ec2go': only_ec2go, 'only_in_obo': only_obo,
           'same_ec_different_go': [{'ec': ec, 'ec2go': sorted(e[ec]), 'obo': sorted(o[ec])} for ec in differ],
           'agreement': {'ecs_in_both': len(set(e) & set(o)), 'same_go': len(set(e) & set(o)) - len(differ)}}
    if a.ec_list:
        ecs = [l.strip() for l in open(a.ec_list) if l.strip()]
        rep['ec_list'] = {'path': a.ec_list, 'n': len(ecs), 'in_ec2go': sum(ec in e for ec in ecs), 'in_obo': sum(ec in o for ec in ecs),
                          'missing_ec2go': [ec for ec in ecs if ec not in e], 'missing_obo': [ec for ec in ecs if ec not in o]}
    Path(a.out).write_text(json.dumps(rep, indent=1) + '\n')
    print(json.dumps({k: v for k, v in rep.items() if k in ('agreement', 'ec_list')}, indent=1), 'only_in_ec2go', len(only_ec2go), 'only_in_obo', len(only_obo), 'differ', len(differ))


if __name__ == '__main__':
    main()
