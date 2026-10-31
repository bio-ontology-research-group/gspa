"""Independent retained-output verification; no Java, native tools or solvers needed."""
from pathlib import Path
import csv,gzip,hashlib,io,json,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent

def sha(data):return hashlib.sha256(data).hexdigest()

def verify():
    for line in (ROOT/'archive.sha256').read_text().splitlines():
        digest,name=line.split('  ',1);assert sha((ROOT/name).read_bytes())==digest,name
    bindings=json.loads((ROOT/'retained-bindings.json').read_text())
    def data(name):
        item=bindings[name];value=(ROOT/item['retained']).read_bytes()
        if item['retained'].endswith('.gz') and not name.endswith('.gz'):value=gzip.decompress(value)
        assert sha(value)==item['sha256'],name
        return value
    def text(name):return data(name).decode()
    def obj(name):return json.loads(text(name))
    for name in bindings:data(name)
    for version in (1,2,3,4):
        manifest=obj(f'launch/repair-source-v{version}.json')
        prefix=f'source-v{version}/'
        for name in bindings:
            if name.startswith(prefix):assert bindings[name]['sha256']==manifest[name[len(prefix):]],name
    def gate(job,phase):
        suites=[ET.fromstring(data(n)) for n in bindings if n.startswith(f'jobs/{job}/{phase}/') and '/TEST-' in n and n.endswith('.xml')]
        assert suites
        return {k:sum(int(s.attrib[k]) for s in suites) for k in ('tests','failures','errors','skipped')}
    assert text('jobs/10513/exit-code.txt').strip()=='23'
    assert text('jobs/10514/exit-code.txt').strip()=='1'
    before=gate('10517','before-tests');assert before==dict(tests=42,failures=23,errors=0,skipped=0),before
    mid=gate('10517','final-tests');assert mid==dict(tests=1548,failures=0,errors=0,skipped=6),mid
    final=gate('10520','final-tests');assert final['tests']>1548 and final['failures']==final['errors']==0 and final['skipped']==6,final
    focused=gate('10520','focused-tests');assert focused['tests']>53 and focused['failures']==focused['errors']==focused['skipped']==0,focused
    for job in ('10515','10517','10520','10521'):assert text(f'jobs/{job}/exit-code.txt').strip()=='0'
    assert 'BUILD SUCCESSFUL' in text('jobs/10520/full.log')
    report=obj('jobs/10521/probe/verification.json');assert report['status']=='verified'
    prefix='jobs/10521/probe/'
    for name,digest in report['output_hashes'].items():assert sha(data(prefix+name))==digest,name
    runs=obj(prefix+'runs.json');assert len(runs)==report['runs']==23
    assert all(r['exit_code']==0 and not r['timed_out'] for r in runs)
    def proteins(name):
        result={}
        for block in text(name).split('>'):
            if block:
                lines=block.splitlines();key=lines[0].split()[0];assert key not in result
                result[key]=''.join(lines[1:])
        return result
    def reference(job,mode):
        base=f'jobs/{job}/observed/genome-{mode}'
        rows=[l.split('\t') for l in text(base+'.gff').splitlines() if l and not l.startswith('#')]
        gff={}
        for row in rows:
            if row[2]!='CDS':continue
            attrs=dict(p.split('=',1) for p in row[8].split(';') if '=' in p)
            assert attrs['ID'] not in gff;gff[attrs['ID']]=row
        result={}
        for block in text(base+'.faa').split('>'):
            if not block:continue
            lines=block.splitlines();h=lines[0].split(' # ')
            attrs=dict(p.split('=',1) for p in h[4].split(';') if '=' in p)
            row=gff.get(h[0],gff.get(attrs['ID']));assert row is not None
            key=(row[0],int(row[3]),int(row[4]),row[6])
            assert (key[1],key[2],key[3])==(int(h[1]),int(h[2]),'+' if h[3]=='1' else '-')
            assert key not in result;result[key]=''.join(lines[1:])
        assert len(gff)==len(result)
        return result
    def output(name,input_name='genome'):
        folder=prefix+name+'/'+input_name
        seq=proteins(folder+'.faa')
        rows=list(csv.DictReader(io.StringIO(text(folder+'_proteins.tsv')),delimiter='\t'))
        assert set(seq)=={r['protein_id'] for r in rows}
        actual={(r['contig'],int(r['start']),int(r['end']),r['strand']):seq[r['protein_id']] for r in rows}
        assert len(actual)==len(rows)
        statuses=obj(folder+'_predictor_status.json')['predictors'];assert len(statuses)==1
        return actual,statuses[0]
    sizes={}
    for tool,job,expected in [('pyrodigal','10507',{'single':4359,'meta':4291}),('prodigal','10515',{'single':4359,'meta':4286})]:
        for mode in ('single','meta'):
            truth=reference(job,mode);assert len(truth)==expected[mode]
            for version in ('old','new'):
                name=version+'-'+tool+'-'+mode
                actual,status=output(name);assert actual==truth
                assert status['outcome']=='SUCCESS'
                if version=='new':assert status['exit_code']==0
            for suffix in ('.faa','_proteins.tsv','.gff3','_annotations.tsv'):
                assert data(prefix+'old-'+tool+'-'+mode+'/genome'+suffix)==data(prefix+'new-'+tool+'-'+mode+'/genome'+suffix)
            sizes[tool+'-'+mode]=len(truth)
        native=obj(f'jobs/{job}/observed/observations.json')
        for mode in ('single','meta'):
            code=next(r['exit_code'] for r in native['cases'] if r['name']=='short-'+mode)
            for version in ('old','new'):
                actual,status=output(version+'-'+tool+'-short-'+mode,'short');assert not actual
                assert status['outcome']==('SUCCESS_ZERO' if version=='new' and mode=='meta' else 'FAILED')
                if version=='new':assert status['exit_code']==code
    reordered,status=output('reordered');assert reordered==reference('10507','single') and status['outcome']=='SUCCESS'
    missing,status=output('missing-protein');assert not missing and status['outcome']=='FAILED'
    assert '\tCDS\t' not in text(prefix+'missing-protein/genome.gff3').split('##FASTA')[0]
    unavailable,status=output('unavailable','short');assert not unavailable and status['outcome']=='UNAVAILABLE'
    for mode,count in [('single',294),('meta',291)]:
        actual,status=output('numeric-'+mode,'numeric');assert len(actual)==count and status['outcome']=='SUCCESS'
        expected={}
        for block in text(f'jobs/10510/observed/numeric-{mode}.faa').split('>'):
            if block:
                lines=block.splitlines();h=lines[0].split(' # ')
                expected[(h[0].rsplit('_',1)[0],int(h[1]),int(h[2]),'+' if h[3]=='1' else '-')]=''.join(lines[1:])
        assert actual==expected
    assert data(prefix+'import-old.jsonl')==data(prefix+'import-new.jsonl')
    claims=[json.loads(l) for l in text(prefix+'import-new.jsonl').splitlines()]
    assert len(claims)==report['imported_unscored_memberships']>0
    assert all(r['metadata']['score_available'] is False for r in claims)
    # Bind every retained runtime input; explicitly enumerate excluded executables/wheels.
    mappings={'/home/leechuck/tmp/gene-probe/observed/':'jobs/10507/observed/',
              '/home/leechuck/tmp/gene-probe/numeric/':'jobs/10510/observed/',
              '/home/leechuck/tmp/gene-probe/prodigal-observed/':'jobs/10515/observed/',
              '/home/leechuck/tmp/gene-probe/results/':prefix,
              '/home/leechuck/Public/software/gspa/':'source-v4/'}
    omitted={};checked=0
    for path,digest in report['bindings'].items():
        if path.endswith('/native.json'):
            assert sha(data('native/antismash.json'))==digest;checked+=1;continue
        mapped=next((dest+path[len(src):] for src,dest in mappings.items() if path.startswith(src)),None)
        if mapped:
            assert sha(data(mapped))==digest,mapped;checked+=1
        else:omitted[path]=digest
    assert len(omitted)==7,omitted # Two JARs, Java, host Python, Prodigal, two wheels.
    assert report['bindings']['/home/leechuck/tmp/gene-probe/after.jar']==text('jobs/10520/jar.sha256').split()[0]
    assert report['bindings']['/home/leechuck/tmp/gene-probe/before.jar']=='069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303'
    return dict(verified=True,retained_files=len(bindings),before=before,focused=focused,full=final,
                actual_cli_runs=len(runs),native_proteins=sizes,unscored_memberships=len(claims),
                retained_runtime_bindings=checked,omitted_runtime_bindings=omitted,
                limitation='Contract and parity verification, not biological accuracy or whole-item acceptance.')

if __name__=='__main__':print(json.dumps(verify(),indent=2))
