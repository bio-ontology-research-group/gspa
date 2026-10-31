#!/usr/bin/env python3
"""Independent native-to-typed comparison and reference EC agreement; no model fitting."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(job):
    assert (job/'exit-code.txt').read_text().strip()=='0', 'native job did not complete'
    before=(job/'native-inputs-before.sha256').read_text().splitlines()
    after=(job/'native-inputs-after.log').read_text().splitlines()
    assert after==[line.split('  ',1)[1]+': OK' for line in before]
    assert json.loads((job/'runtime.json').read_text())['cuda_available'] is True
    result={}
    for dataset,n_expected in [('new',392),('price',149)]:
        assert (job/(dataset+'.exit-code.txt')).read_text().strip()=='0'
        with (job/'truth'/f'{dataset}.csv').open() as stream:
            records=list(csv.DictReader(stream,delimiter='\t'))
        truth={r['Entry']:{'EC:'+ec.strip() for ec in r['EC number'].split(';') if ec.strip()} for r in records}
        assert len(records)==len(truth)==n_expected
        assert all(truth.values())
        sequences={};pid=None
        for line in (job/'inputs'/f'{dataset}.fasta').read_text().splitlines():
            if line.startswith('>'):
                pid=line[1:].split()[0];assert pid not in sequences;sequences[pid]=''
            elif line.strip():
                assert pid is not None;sequences[pid]+=line.strip()
        assert sequences=={r['Entry']:r['Sequence'] for r in records},'input differs from primary reference sequences'
        folder=job/'predictions'/dataset
        native=folder/(dataset+'.clean.native.csv')
        tsv=folder/(dataset+'.clean.tsv')
        receipt=json.loads((folder/(dataset+'.clean.receipt.json')).read_text())
        assert receipt['input_fasta_sha256']==sha(job/'inputs'/f'{dataset}.fasta')
        assert receipt['native_output_sha256']==sha(native)
        assert receipt['output_sha256']==sha(tsv)
        assert receipt['min_score']==0 and receipt['filtered_calls']==0
        assert receipt['score_kind']=='gmm_confidence' and receipt['score_available'] is True
        assert set(receipt['input_protein_ids'])==set(truth)
        mapping=receipt['native_to_input_ids']
        assert len(mapping)==n_expected and set(mapping.values())==set(truth)
        calls={};proteins=set()
        with native.open() as stream:
            for row in csv.reader(stream):
                assert row[0] not in proteins and len(row)>1
                proteins.add(row[0]);pid=mapping[row[0]]
                for token in row[1:]:
                    term,value=token.split('/')
                    score=float(value)
                    assert math.isfinite(score) and 0<=score<=1 and (pid,term) not in calls
                    calls[pid,term]=(score,value,row[0])
        assert proteins==set(mapping)
        typed={}
        with tsv.open() as stream:
            header=stream.readline().rstrip('\n').split('\t')
            assert header==['protein_id','term','score','annotation_type','metadata_json']
            for line in stream:
                pid,term,score,kind,metadata=line.rstrip('\n').split('\t')
                metadata=json.loads(metadata);key=(pid,term)
                assert key in calls and key not in typed and kind=='EC'
                value,text,native_id=calls[key]
                assert float(score)==metadata['native_score']==value
                assert metadata['native_score_text']==text and metadata['native_protein_id']==native_id
                assert metadata['score_available'] is True and metadata['score_calibrated'] is False
                assert metadata['native_score_direction']=='higher' and metadata['score_kind']=='gmm_confidence'
                typed[key]=value
        assert set(typed)==set(calls)
        assert receipt['native_proteins']==receipt['input_proteins']==n_expected
        assert receipt['native_selected_calls']==receipt['emitted_calls']==len(calls)
        predicted={pid:set() for pid in truth}
        for pid,term in calls:predicted[pid].add(term)
        hits=sum(len(truth[pid]&predicted[pid]) for pid in truth)
        labels=sum(map(len,truth.values()));selected=len(calls)
        precision=hits/selected;recall=hits/labels
        result[dataset]={'proteins':n_expected,'native_selected_calls':selected,'reference_labels':labels,
            'reference_matches':hits,'predictions_absent_from_reference':selected-hits,
            'unrecovered_reference_labels':labels-hits,'reference_precision':precision,'reference_recall':recall,
            'reference_micro_f1':2*precision*recall/(precision+recall) if precision+recall else 0,
            'exact_reference_set_proteins':sum(predicted[pid]==truth[pid] for pid in truth),
            'any_reference_match_proteins':sum(bool(predicted[pid]&truth[pid]) for pid in truth),
            'zero_native_confidences':sum(v[0]==0 for v in calls.values()),
            'provisional_ec_calls':sum('.n' in term for pid,term in calls),
            'sha256':{'native':sha(native),'typed':sha(tsv),'input':sha(job/'inputs'/f'{dataset}.fasta'),
                      'reference':sha(job/'truth'/f'{dataset}.csv')}}
    return {'status':'complete_native_and_typed_agreement','datasets':result,'proteins':541,
        'scope':'CLEAN maxsep selected calls at zero cutoff; reference agreement, not independent functional accuracy or GSPA improvement.',
        'limitations':['Reference absence is not an experimentally established negative.',
                       'Training/pretraining overlap not assessed; installed model files recorded by SHA256.',
                       'No complete EC distance matrix, calibrated probability claim, context comparison or function-recovery claim.']}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('job',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=audit(args.job)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,sort_keys=True))
