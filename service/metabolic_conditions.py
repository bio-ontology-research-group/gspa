"""Test nutrient dependence on the exported model, without reconstructing it again."""
import argparse
import json
from pathlib import Path
import libsbml
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

AMINO_ACIDS = dict(zip(
    ['00035','00051','00132','00041','00084','00053','00023','00033','00119','00322',
     '00107','00039','00060','00066','00129','00054','00161','00065','00069','00156'],
    ['alanine','arginine','asparagine','aspartate','cysteine','glutamine','glutamate','glycine','histidine','isoleucine',
     'leucine','lysine','methionine','phenylalanine','proline','serine','threonine','tryptophan','tyrosine','valine']))


def evaluate(path):
    document=libsbml.readSBML(str(path))
    if document.getNumErrors(): raise ValueError(document.getErrorLog().toString())
    model=document.getModel()
    mets={s.id:i for i,s in enumerate(model.getListOfSpecies())}
    reactions=list(model.getListOfReactions()); indices={r.id:i for i,r in enumerate(reactions)}
    ri=[];ci=[];values=[];bounds=[]
    for col,r in enumerate(reactions):
        fbc=r.getPlugin('fbc')
        bounds.append((model.getParameter(fbc.getLowerFluxBound()).value,model.getParameter(fbc.getUpperFluxBound()).value))
        for sign,refs in ((-1,r.getListOfReactants()),(1,r.getListOfProducts())):
            for ref in refs:ri.append(mets[ref.species]);ci.append(col);values.append(sign*ref.stoichiometry)
    matrix=coo_matrix((values,(ri,ci)),shape=(len(mets),len(reactions))).tocsc()
    objective=model.getPlugin('fbc').getActiveObjective()
    cost=np.zeros(len(reactions))
    for term in objective.getListOfFluxObjectives():
        col=indices[term.reaction];cost[col]=-term.coefficient
        bounds[col]=(0,bounds[col][1])  # remove the reconstruction target to permit no-growth solutions
    conditions=[('Supplemented reference medium',[])]
    conditions += [('Without '+name+' uptake',['EX_cpd'+cpd+'_e']) for cpd,name in AMINO_ACIDS.items()]
    conditions += [('Without the 20 listed amino-acid uptakes',['EX_cpd'+cpd+'_e' for cpd in AMINO_ACIDS])]
    results=[]
    for label,closed in conditions:
        changed=list(bounds)
        for rid in closed:
            if rid not in indices:continue  # absent exchange already cannot import this metabolite
            col=indices[rid];coeff=matrix[:,col].data
            if len(coeff)!=1: raise ValueError('Expected a one-metabolite exchange: '+rid)
            lo,hi=changed[col]
            changed[col]=(lo,min(0,hi)) if coeff[0]>0 else (max(0,lo),hi)
        result=linprog(cost,A_eq=matrix,b_eq=np.zeros(len(mets)),bounds=changed,method='highs')
        if not result.success and result.status!=2:raise RuntimeError('Nutrient test failed: '+result.message)
        flux=max(0,float(-result.fun)) if result.success else 0
        results.append({'condition':label,'closed_exchanges':closed,'maximum_biomass_flux':flux,
                        'status':'feasible biomass flux' if flux>1e-7 else 'no feasible biomass flux'})
    return {'method':'Maximize biomass on the exported network. Each test closes only the indicated uptake exchanges; other medium bounds remain unchanged. No additional gap filling.',
            'interpretation':'Conditional model predictions, not observed growth or validated auxotrophies. Model omissions can cause apparent dependencies.',
            'results':results}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.write_text(json.dumps(evaluate(args.model),indent=2))
