#!/usr/bin/env python3
"""Recheck retained native calls, complete operon partitions and controlled route scores."""
import csv
import gzip
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def content(path):
    return gzip.open(path,'rb').read() if path.suffix=='.gz' else path.read_bytes()


def protein_sequences(path):
    rows={};key=None
    for line in content(path).decode().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0];assert key not in rows;rows[key]=''
        elif line:
            assert key is not None;rows[key]+=line.strip()
    return rows


def coordinates(path):
    rows={}
    for f in csv.reader(io.StringIO(content(path).decode()),delimiter='\t'):
        if not f or f[0].startswith('#'):continue
        assert len(f)==9 and f[2]=='CDS'
        attrs=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
        key=attrs['ID'];assert key not in rows
        rows[key]=(f[0],int(f[3]),int(f[4]),f[6],f[7])
    return rows


def score_rows(path):
    return {(r['protein_id'],r['function_id']):r for r in csv.DictReader(path.open(),delimiter='\t')}


def main():
    entries=(ROOT/'archive.sha256').read_text().splitlines()
    for line in entries:
        digest,name=line.split('  ',1)
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    for version in ('v1','v2'):
        for name,digest in json.loads((ROOT/f'route-source-{version}.json').read_text()).items():
            assert hashlib.sha256((ROOT/f'source-{version}'/name).read_bytes()).hexdigest()==digest,name
    bindings=json.loads((ROOT/'bindings-10490.json').read_text())
    for name,item in bindings.items():
        if 'path' in item:assert hashlib.sha256(content(ROOT/item['path'])).hexdigest()==item['sha256'],name
    assert sum('path' in v for v in bindings.values())==44
    assert sum('omitted' in v for v in bindings.values())==5
    cli=ROOT/'jobs/10490/cli';observed=ROOT/'observed/10484'
    assert (ROOT/'jobs/10490/exit-code.txt').read_text().strip()=='0'
    report=json.loads((cli/'report.json').read_text())
    assert len(report['commands'])==9 and all(c['exit_code']==0 for c in report['commands'])
    assert report['inputs_sha256']=={name:v['sha256'] for name,v in bindings.items()}
    assert bindings['gspa-cli/build/libs/gspa-3.0.0.jar']['sha256']=='069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303'
    assert hashlib.sha256(content(observed/'genome.fna.gz')).hexdigest()=='76a0087d8d90585696bd21d4c01b0e1f6699c72cce2a9a4682898c0255154ed5'
    result={}
    for mode,nprot,nops,nmembers in [('single',4359,819,3867),('meta',4291,812,3789)]:
        sample=cli/(mode+'-output')/mode
        proteins=protein_sequences(sample/'gene_calling'/(mode+'.faa.gz'))
        loci=coordinates(sample/'gene_calling'/(mode+'.gff.gz'))
        assert len(proteins)==nprot and proteins.keys()==loci.keys()
        assert proteins==protein_sequences(observed/(mode+'.faa.gz'))
        assert loci==coordinates(observed/(mode+'.gff.gz'))
        opsfile=sample/(mode+'_operons_for_integrate.tsv')
        groups=[line.split('\t') for line in opsfile.read_text().splitlines()]
        assert opsfile.read_bytes()==(observed/(mode+'-operons-before/operons_for_integrate.tsv')).read_bytes()
        # Independently construct connected components of all qualifying adjacency edges.
        parent={pid:pid for pid in loci};by_contig=defaultdict(list)
        def find(pid):
            while parent[pid]!=pid:
                parent[pid]=parent[parent[pid]];pid=parent[pid]
            return pid
        for pid,loc in loci.items():by_contig[loc[0]].append(pid)
        for members in by_contig.values():
            members.sort(key=lambda pid:loci[pid][1])
            for a,b in zip(members,members[1:]):
                if loci[a][3]==loci[b][3] and loci[b][1]-loci[a][2]-1<=300:
                    parent[find(b)]=find(a)
        components=defaultdict(set)
        for pid in loci:components[find(pid)].add(pid)
        expected={frozenset(c) for c in components.values() if len(c)>=2}
        assert len(groups)==nops and sum(map(len,groups))==nmembers
        assert len({frozenset(g) for g in groups})==len(groups)
        assert {frozenset(g) for g in groups}==expected
        assert len([pid for g in groups for pid in g])==len({pid for g in groups for pid in g})
        integrated=sample/(mode+'_integrated.tsv')
        assert integrated.read_bytes()==(cli/(mode+'-direct.tsv')).read_bytes()
        live,base,bad=map(score_rows,[integrated,cli/(mode+'-baseline.tsv'),cli/(mode+'-numeric.tsv')])
        assert live.keys()==base.keys()==bad.keys() and len(live)==4
        summary=report['summary'][mode];candidate=(summary['candidate'],'GO:0004332')
        assert 0.18<float(base[candidate]['posterior_prob'])<0.5<float(live[candidate]['posterior_prob'])
        assert 'genomic_context' in live[candidate]['priors_fired']
        for key in live:
            assert float(base[key]['posterior_prob'])==float(bad[key]['posterior_prob'])
            if key!=candidate:
                assert all(live[key][field]==base[key][field] for field in ('posterior_prob','likelihood_logodds','final_logodds','n_supporting','priors_fired'))
        assert (summary['proteins'],summary['operons'],summary['operon_members'])==(nprot,nops,nmembers)
        result[mode]=dict(proteins=nprot,complete_operon_groups=nops,members=nmembers)
    print(json.dumps(dict(files=len(entries),retained_bindings=44,omitted_runtime_bindings=5,verified=result),indent=2))


if __name__=='__main__':main()
