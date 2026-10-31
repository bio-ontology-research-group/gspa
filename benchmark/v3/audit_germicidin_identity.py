"""Compare measured germicidin identities with the pinned model, without relabelling.

A matching name or molecular formula is not proof of structural identity. This
source audit produces no function labels and never changes model chemistry.
"""
import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

SBML='http://www.sbml.org/sbml/level3/version1/core'
FBC='http://www.sbml.org/sbml/level3/version1/fbc/version2'
NS={'s':SBML,'f':FBC}
REACTIONS=('R_GERMA','R_GERMB','R_GERMC','R_GERMD')
PINNED_ARTICLE='65159aa2bf3fa6717b47d7c78b47825d78ee9b0edb8371fd3c3c1fc362a699f1'
PINNED_MODEL='2d2798faf43bebba0f9b4fffe147d68c9050276eaaf28bc483f90d99c5b107b6'
HEADER=['Secondary metabolite','Medium (carbon source)','Molecular formula',
        'Theoretical mass (m/z)','Acquired mass (m/z)','Mass error (ppm)']

class Tables(HTMLParser):
    def __init__(self):
        super().__init__();self.active=False;self.tables={};self.cell=None;self.table_id=None
        self.paragraph=None;self.paragraphs=[]
    def handle_starttag(self,tag,attrs):
        attr=dict(attrs)
        if re.fullmatch(r'T[1-3]',attr.get('id','')):self.table_id=attr['id']
        if tag=='p':self.paragraph=[]
        if self.active and any(attr.get(k,'1')!='1' for k in ('colspan','rowspan')):
            raise ValueError('spanning source cell unsupported')
        if tag=='table':
            if self.active:raise ValueError('nested source table')
            if self.table_id is None or self.table_id in self.tables:raise ValueError('missing or duplicate table ID')
            self.active=True;self.current=[]
        if self.active and tag=='tr':self.row=[]
        if self.active and tag in ('td','th'):self.cell=[]
    def handle_data(self,data):
        if self.active and self.cell is not None:self.cell.append(data)
        if self.paragraph is not None:self.paragraph.append(data)
    def handle_endtag(self,tag):
        if tag=='p' and self.paragraph is not None:
            self.paragraphs.append(''.join(self.paragraph).strip());self.paragraph=None
        if self.active and tag in ('td','th'):
            self.row.append(''.join(self.cell).strip());self.cell=None
        if self.active and tag=='tr':self.current.append(self.row)
        if self.active and tag=='table':self.tables[self.table_id]=self.current;self.active=False

def composition(formula):
    pieces=re.findall(r'([A-Z][a-z]?)([0-9]*)',formula)
    if ''.join(e+n for e,n in pieces)!=formula or not pieces:raise ValueError('unsupported formula')
    result={}
    for element,count in pieces:
        n=int(count or '1')
        if n<=0:raise ValueError('nonpositive atom count')
        result[element]=result.get(element,0)+n
    return result

def model_products(path):
    root=ET.parse(path)
    metabolites={s.attrib['id']:s for s in root.findall('.//s:listOfSpecies/s:species',NS)}
    reactions={r.attrib['id']:r for r in root.findall('.//s:listOfReactions/s:reaction',NS)}
    out=[]
    for rid in REACTIONS:
        reaction=reactions[rid]
        products=[x.attrib['species'] for x in reaction.findall('s:listOfProducts/s:speciesReference',NS)
                  if x.attrib['species'].startswith('M_germicidin')]
        if len(products)!=1:raise ValueError('unexpected germicidin product population')
        metabolite=metabolites[products[0]];formula=metabolite.attrib['{'+FBC+'}chemicalFormula']
        atoms=composition(formula)
        if not set(atoms)<=set('CHO'):raise ValueError('unexpected product chemistry')
        out.append(dict(reaction_id=rid,metabolite_id=products[0],name=metabolite.attrib['name'],
            formula=formula,composition=atoms,charge=int(metabolite.attrib['{'+FBC+'}charge']),
            metabolite_annotations=[x.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource']
                for x in metabolite.findall('.//{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li')],
            nominal_mass_assumption='Derived neutral CHO formula plus one proton: nominal mass, not exact m/z.',
            nominal_protonated_mass=sum(atoms.get(e,0)*w for e,w in [('C',12),('H',1),('O',16)])+1,
            reactants=[dict(metabolite_id=x.attrib['species'],name=metabolites[x.attrib['species']].attrib['name'],
                            stoichiometry=x.attrib['stoichiometry'])
                       for x in reaction.findall('s:listOfReactants/s:speciesReference',NS)],
            model_gene_rules=[x.attrib['{'+FBC+'}geneProduct'] for x in reaction.findall('.//f:geneProductRef',NS)],
            annotation_links=[x.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource']
                              for x in reaction.findall('.//{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li')]))
    return out

def audit(article,model):
    for path,expected in [(article,PINNED_ARTICLE),(model,PINNED_MODEL)]:
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('pinned input changed: '+str(path))
    html=article.read_text(encoding='utf-8')
    parser=Tables();parser.feed(html)
    if set(parser.tables)!={'T1','T2','T3'} or any(parser.tables[t][0]!=HEADER for t in ['T2','T3']):
        raise ValueError('source table layout changed')
    if parser.tables['T1'][0]!=['Secondary metabolites','Biosynthetic (SCO) genes expressed during germination']:
        raise ValueError('source biosynthesis header changed')
    if any(len(row)!=2 for row in parser.tables['T1']):raise ValueError('source biosynthesis width changed')
    tables={key[1:]:rows for key,rows in parser.tables.items()}
    def paragraph(phrase):
        # The publisher repeats identical captions in its navigation and body.
        matches=list(dict.fromkeys(p for p in parser.paragraphs if phrase in p))
        if len(matches)!=1:raise ValueError('missing or duplicate primary context: '+phrase)
        return matches[0]
    context={
        '2':dict(assay_context='stationary_phase',caption=paragraph('Secondary metabolites detected in the stationary phase of growth.'),
                 standard_confirmation='not separately established by this audit'),
        '3':dict(assay_context='germination',caption=paragraph('Secondary metabolites produced during germination.'),
                 identity_evidence=paragraph('confirmed by comparing the actual retention time with the original standard'))}
    figure=paragraph('Acquired m/z 197.1188 corresponds to')
    caption_pair=re.search(r'Acquired m/z (197\.[0-9]+).*?the theoretical m/z (197\.[0-9]+)',figure)
    if caption_pair is None:raise ValueError('primary figure mass description changed')
    measured=[]
    for number,rows in [('2',parser.tables['T2']),('3',parser.tables['T3'])]:
        for index,row in enumerate(rows[1:],1):
            if len(row)!=6:raise ValueError('source row width changed')
            if row[0].startswith('Germicidin'):
                measured.append(dict(table=number,data_row=index,**dict(zip(HEADER,row))))
    if [(r['table'],r['Secondary metabolite']) for r in measured]!=[('2','Germicidin A'),('2','Germicidin B'),('3','Germicidin A')]:
        raise ValueError('measured congener population changed')
    products=model_products(model);comparisons=[]
    for row in measured:
        atoms=composition(row['Molecular formula'])
        same_name=[p for p in products if p['name'].casefold()==row['Secondary metabolite'].casefold()]
        same_formula=[p['metabolite_id'] for p in products if p['composition']==atoms and p['charge']==0]
        disposition=('no_formula_compatible_candidate' if not same_formula else
            'name_formula_conflict' if same_name and any(p['composition']!=atoms for p in same_name) else
            'formula_ambiguous' if len(same_formula)>1 else 'formula_compatible_structure_unverified')
        comparisons.append(dict(observation=row,source_context=context[row['table']],name_matches=[dict(metabolite_id=p['metabolite_id'],
            model_formula=p['formula'],formula_agrees=p['composition']==atoms,
            model_nominal_protonated_mass=p['nominal_protonated_mass']) for p in same_name],
            formula_compatible_candidates=same_formula,accepted_identity=None,
            disposition=disposition))
    return dict(schema='gspa.germicidin-identity-audit/1',
        inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [article,model]},
        primary_tables=tables,model_products=products,comparisons=comparisons,
        source_consistency=dict(figure2_caption=figure,figure2_panel_b_acquired_mz=caption_pair.group(1),
            figure2_panel_b_theoretical_mz=caption_pair.group(2),
            table3_acquired_mz=measured[2]['Acquired mass (m/z)'],table3_theoretical_mz=measured[2]['Theoretical mass (m/z)'],
            table3_caption_labels_disagree=(caption_pair.group(1)!=measured[2]['Acquired mass (m/z)'] or
                caption_pair.group(2)!=measured[2]['Theoretical mass (m/z)']),
            interpretation='Reported values retained, not repaired or treated as structural identification; no exact-mass or ppm validation.'),
        scope='Source identity audit only; no empirical labels, model relabelling, or accuracy calculation.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('article','model','out'):p.add_argument('--'+key,type=Path,required=True)
    args=p.parse_args();report=audit(args.article,args.model)
    serialized=json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n'
    with args.out.open('x',encoding='utf-8') as f:f.write(serialized)
