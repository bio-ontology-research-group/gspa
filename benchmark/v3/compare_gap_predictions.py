"""Measure full-workflow gap effects without calling them correct or harmful.

Gaps activate the existing gap-filling prior as well as RLC. Identity and sequence
likelihood must remain unchanged; finite posterior/log-odds changes are retained
for later independent evaluation, including decreases and threshold crossings.
"""
from collections import Counter
import csv
import math
from pathlib import Path


def read(path):
    result={}
    with Path(path).open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        for row in reader:
            key=tuple(row[k] for k in ('protein_id','type','function_id'))
            if key in result:raise ValueError('duplicate prediction identity')
            if None in row or None in row.values():raise ValueError('malformed prediction row')
            result[key]=row
    return result


def compare(before, after, differences):
    a,b=read(before),read(after)
    if a.keys()!=b.keys():raise ValueError('functional prediction identity changed')
    counts=Counter();delta_rows=[];maximum=0.
    for key,old in a.items():
        new=b[key]
        if old.keys()!=new.keys():raise ValueError('prediction columns changed')
        for field in ('go_aspect','likelihood_logodds','n_supporting'):
            if old[field]!=new[field]:raise ValueError('sequence evidence changed: '+field)
        for row in (old,new):
            prob=float(row['posterior_prob']);logodds=float(row['final_logodds']);likelihood=float(row['likelihood_logodds'])
            if not all(map(math.isfinite,(prob,logodds,likelihood))) or not 0<=prob<=1:
                raise ValueError('nonfinite or invalid posterior')
            expected=1/(1+math.exp(-logodds)) if logodds>=0 else math.exp(logodds)/(1+math.exp(logodds))
            if abs(expected-prob)>1e-12:raise ValueError('posterior/log-odds mismatch')
        x,y=float(old['posterior_prob']),float(new['posterior_prob']);delta=y-x
        counts['predictions']+=1
        if delta:
            counts['changed']+=1;counts['increased' if delta>0 else 'decreased']+=1
            counts['crossed_above_0.5']+=int(x<.5<=y)
            counts['crossed_below_0.5']+=int(y<.5<=x)
            maximum=max(maximum,abs(delta))
            delta_rows.append([*key,x,y,delta,old['priors_fired'],new['priors_fired']])
    with Path(differences).open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['protein_id','type','function_id','without_gaps','with_gaps','difference','priors_without_gaps','priors_with_gaps'])
        w.writerows(delta_rows)
    return dict(counts=dict(counts),max_absolute_change=maximum,
                sequence_likelihood_unchanged=True,threshold_interpretation='Descriptive crossings at 0.5; not calibrated or biologically verified decisions.')
