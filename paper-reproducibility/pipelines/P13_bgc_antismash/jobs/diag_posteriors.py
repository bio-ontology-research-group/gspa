#!/usr/bin/env python3
"""Diagnostic: posterior distribution of class-term claims on core genes, per configuration (P13)."""
import sys, json, os, importlib.util, statistics
from collections import defaultdict, Counter
O0 = os.environ.get("OUT", "/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19")
spec=importlib.util.spec_from_file_location("b",O0 + "/code/bgc_benchmark.py"); b=importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
O=O0
n,ns,par,obs=b.read_obo(O+"/ref/go.obo")
terms={c:b.CLASS_MAP[c][0] for c in b.CLASS_ORDER if c!="ALL_BGC"}
desc=b.descendants(par,set(terms.values()))
rev=defaultdict(list)
for r,ds in desc.items():
    for d in ds: rev[d].append(r)
gids=[l.split("\t")[0] for l in open(O+"/genomes.tsv")][1:]
post={c:[] for c in ["stage0_noisy_or","operon_only","gspa_A9"]}
src=Counter(); rawmax=[]; changed=0; tot=0
for g in gids:
    genes={}
    for l in open(f"{O}/work/{g}/{g}_genes.tsv").readlines()[1:]:
        f=l.rstrip("\n").split("\t")
        if f[9]=="core": genes[f[1]]=set(f[7].split(","))
    best={}
    for l in open(f"{O}/work/{g}/{g}_claims.jsonl"):
        c=json.loads(l)
        p=c["protein_id"]
        if p not in genes: continue
        for r in rev.get(c["function_id"],()):
            cl=[k for k,v in terms.items() if v==r][0]
            if cl in genes[p]:
                src[c["source"]]+=1
                best[(p,r)]=max(best.get((p,r),0),c["raw_score"])
    vals={}
    for cfg in post:
        d={}
        for l in open(f"{O}/integrated/{cfg}/{g}.tsv").readlines()[1:]:
            f=l.split("\t")
            for r in rev.get(f[2],()):
                k=(f[0],r)
                if k in best: d[k]=max(d.get(k,0),float(f[4]))
        for k in best: post[cfg].append(d.get(k,0.0))
        vals[cfg]=d
    for k in best:
        tot+=1
        if abs(vals["gspa_A9"].get(k,0)-vals["stage0_noisy_or"].get(k,0))>0.01: changed+=1
print("core (gene,class-term) pairs with >=1 claim:",tot,"; changed by >0.01 under A9:",changed)
print("claims by source (class term or descendant, core genes):",dict(src))
for cfg,v in post.items():
    v=sorted(v); q=lambda x: v[int(x*(len(v)-1))]
    print(cfg,"posterior quantiles 10/25/50/75/90:",[round(q(x),3) for x in (.1,.25,.5,.75,.9)],"frac>=0.5:",round(sum(x>=0.5 for x in v)/len(v),3),"frac>=0.3:",round(sum(x>=0.3 for x in v)/len(v),3))
