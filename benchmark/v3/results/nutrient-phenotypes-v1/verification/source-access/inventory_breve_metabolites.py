from pathlib import Path
import collections,json,re,hashlib
import openpyxl
b=Path.home()/'.cache/gspa-v3-data/breve-functional-genomics-v1'
p=b/'raw/aromatic_lactic_acid/230117_Metabolomics/data/KCHU001_PreferredMetabolitesList.xlsx'
w=openpyxl.load_workbook(p,read_only=True,data_only=True)
rows=[list(r) for r in w['Sheet2'].iter_rows(min_row=9,values_only=True) if any(v is not None for v in r)]
counts=collections.Counter(str(r[0]) for r in rows)
m=b/'raw/essentialgaps/data/name.conversion.xlsx'
mapping=collections.defaultdict(list)
for name,rs,old in openpyxl.load_workbook(m,read_only=True,data_only=True).active.iter_rows(min_row=2,values_only=True):
 if rs:mapping[rs].append(dict(symbol=name,old=old))
resolved=[]
for sample,n in sorted(counts.items()):
 match=re.fullmatch(r'(RS\d+)_P(\d+)_([A-H]\d+)',sample)
 resolved.append(dict(sample=sample,rows=n,canonical='BBR_'+match[1] if match else None,mapping=mapping.get('BBR_'+match[1],[]) if match else []))
report=dict(status='inventory_only_no_labels',source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),mapping_sha256=hashlib.sha256(m.read_bytes()).hexdigest(),nonempty_rows=len(rows),sample_groups=len(counts),samples=resolved,
 author_processing='Mean media blank subtraction per metabolite, floor to 1, log2 transform, mean by strain, subtract pooled-mutant control mean; not absolute concentrations or WT normalization.',
 author_uncertainty='Code combines strain and pool SD in quadrature then divides by sqrt(strain n). Pool n is not independently propagated; do not reuse as a validated independent-sample standard error.',
 caution=['Author code labels within-strain repeats technical; biological independence must be audited before uncertainty claims.','Four primary products have different MSI support; HICA annotation includes MSI 1,4.','Insertion mutants can affect neighbouring genes.','Primary model was curated using fitness data; no independent predictive validity established.'])
(b/'metabolite-inventory.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='samples'},indent=2))
print('focus', [x for x in resolved if x['canonical'] in ['BBR_RS15505','BBR_RS15550'] or x['sample'] in ['pool','media_blank']])
