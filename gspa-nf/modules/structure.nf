/*
 * STRUCTURE_PROVIDER — produce a per-sample directory of PDB models for
 * the predictors that need structures (ScanNet, DeepFRI structure mode).
 *
 * Modes (params.structures_from):
 *   - 'esmfold' : run ESMFold over each protein (5 sec/protein on GPU)
 *   - 'afdb'    : look up UniProt accession → AlphaFold DB model (curl)
 *   - 'none'    : skip; downstream structure consumers will skip with warn
 *
 * Output: directory of <protein_id>.pdb files for the sample.
 */

process STRUCTURE_PROVIDER {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/structures", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}_structures"), emit: structures

    script:
    """
    mkdir -p ${sample_id}_structures
    if [[ "${params.structures_from}" == "esmfold" ]]; then
        python3 -c "
import sys, torch, esm
from pathlib import Path
model = esm.pretrained.esmfold_v1().eval()
if torch.cuda.is_available():
    model = model.cuda()

def fasta(p):
    name=None; chunks=[]
    for line in open(p):
        line=line.rstrip()
        if line.startswith('>'):
            if name: yield name, ''.join(chunks)
            name=line[1:].split()[0]; chunks=[]
        else: chunks.append(line)
    if name: yield name, ''.join(chunks)

out = Path('${sample_id}_structures')
out.mkdir(exist_ok=True)
for pid, seq in fasta('${proteins}'):
    if not seq: continue
    seq = seq[:1024]  # length cap
    with torch.no_grad():
        pdb = model.infer_pdb(seq)
    (out / f'{pid}.pdb').write_text(pdb)
    print(f'  wrote {pid}.pdb', flush=True)
"
    elif [[ "${params.structures_from}" == "afdb" ]]; then
        # Map FASTA headers to UniProt accessions then curl AFDB
        grep '^>' ${proteins} | awk '{print substr(\$1,2)}' | while read pid; do
            # Naive: try as UniProt accession; skip on 404
            url="https://alphafold.ebi.ac.uk/files/AF-\${pid}-F1-model_v4.pdb"
            curl -fsSL -o ${sample_id}_structures/\${pid}.pdb "\$url" 2>/dev/null \
                || echo "  no AFDB hit for \$pid" >&2
        done
    else
        echo "structures_from=none — skipping (downstream structure consumers will skip)"
    fi
    echo "structures: \$(ls ${sample_id}_structures | wc -l) PDBs"
    """
}
