import os,sys
from pathlib import Path

def infer_maxsep(training,name,gmm=None):
 ids=[line[1:].split()[0] for line in Path('data/inputs/'+name+'.fasta').read_text().splitlines() if line.startswith('>')]
 fault=os.environ.get('CLEAN_FAULT')
 if fault=='missing':return
 if fault=='unknown':ids[0]='unknown'
 values=('0.00000001','0.8') if gmm else ('0.01','8.0627')
 text=''.join(pid+',EC:1.1.1.1/'+values[0]+',EC:1.1.1.n11/'+values[1]+'\n' for pid in ids)
 if fault=='malformed':text='malformed\n'
 Path('results/inputs/'+name+'_maxsep.csv').write_text(text)

infer_maxsep('fixture',sys.argv[-1],gmm=None)
