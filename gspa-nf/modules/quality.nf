process MERGE_ANNOTATIONS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'python:3.12-slim'

    input:
    tuple val(sample_id), path(gff), path(proteins), path(diamond_tsv), path(pfam_domtbl), path(interproscan_tsv), path(rrna_gff), path(crispr_gff)

    output:
    tuple val(sample_id), path("${sample_id}_annotations.tsv"), emit: annotations
    tuple val(sample_id), path("${sample_id}_annotated.gff3"),  emit: gff
    tuple val(sample_id), path("${sample_id}_summary.txt"),     emit: summary

    script:
    """
    #!/usr/bin/env python3
    import os, re
    from collections import defaultdict

    annotations = defaultdict(list)

    # Parse DIAMOND results
    diamond_file = "${diamond_tsv}"
    if os.path.exists(diamond_file) and os.path.getsize(diamond_file) > 0:
        with open(diamond_file) as f:
            for line in f:
                fields = line.strip().split('\\t')
                if len(fields) < 7: continue
                qid, sid, pident = fields[0], fields[1], fields[2]
                evalue = fields[6] if len(fields) > 6 else ''
                title = fields[8] if len(fields) > 8 else ''
                acc = sid.split('|')[1] if '|' in sid else sid
                m = re.search(r'(?:sp|tr)\\|\\S+\\|\\S+\\s+(.+?)\\s+OS=', title)
                product = m.group(1) if m else ''
                for go in re.findall(r'GO:\\d{7}', title):
                    annotations[qid].append(('GO', go, pident, 'diamond', 'IEA'))
                for ec in re.findall(r'EC:[\\d\\-]+\\.[\\d\\-]+\\.[\\d\\-]+\\.[\\d\\-]+', title):
                    annotations[qid].append(('EC', ec, pident, 'diamond', 'IEA'))
                if product:
                    annotations[qid].append(('product', product, pident, 'diamond', 'ISS'))

    # Parse HMMER domtblout
    pfam_file = "${pfam_domtbl}"
    if os.path.exists(pfam_file) and os.path.getsize(pfam_file) > 0:
        with open(pfam_file) as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.split()
                if len(fields) < 22: continue
                protein_id = fields[0]
                pfam_acc = re.sub(r'\\.\\d+\$', '', fields[4])
                pfam_name = fields[3]
                dom_score = fields[13]
                annotations[protein_id].append(('Pfam', pfam_acc, dom_score, 'hmmer', pfam_name))

    # Parse barrnap rRNA
    rrna_file = "${rrna_gff}"
    if os.path.exists(rrna_file):
        with open(rrna_file) as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                fields = line.strip().split('\\t')
                if len(fields) >= 9 and fields[2] == 'rRNA':
                    attrs = dict(kv.split('=',1) for kv in fields[8].split(';') if '=' in kv)
                    product = attrs.get('product', attrs.get('Name', 'rRNA'))
                    key = f"{fields[0]}:{fields[3]}-{fields[4]}"
                    annotations[key].append(('rRNA', product, fields[5], 'barrnap', ''))

    # Write merged TSV
    with open("${sample_id}_annotations.tsv", 'w') as out:
        out.write('protein_id\\ttype\\tvalue\\tscore\\tsource\\tevidence\\n')
        for pid in sorted(annotations.keys()):
            for ann_type, value, score, source, evidence in annotations[pid]:
                out.write(f'{pid}\\t{ann_type}\\t{value}\\t{score}\\t{source}\\t{evidence}\\n')

    # Build annotation lookup
    protein_go = defaultdict(set)
    protein_ec = defaultdict(set)
    protein_pfam = defaultdict(set)
    protein_product = {}
    for pid, anns in annotations.items():
        for ann_type, value, score, source, evidence in anns:
            if ann_type == 'GO': protein_go[pid].add(value)
            elif ann_type == 'EC': protein_ec[pid].add(value)
            elif ann_type == 'Pfam': protein_pfam[pid].add(value)
            elif ann_type == 'product' and pid not in protein_product:
                protein_product[pid] = value

    # Create annotated GFF3
    with open("${gff}") as fin, open("${sample_id}_annotated.gff3", 'w') as fout:
        fout.write('##gff-version 3\\n')
        for line in fin:
            if line.startswith('#'):
                if not line.startswith('##gff-version'):
                    fout.write(line)
                continue
            fields = line.strip().split('\\t')
            if len(fields) < 9: continue
            attrs = dict(kv.split('=',1) for kv in fields[8].split(';') if '=' in kv)
            pid = attrs.get('ID', '')
            extra = []
            if pid in protein_go:
                extra.append('Ontology_term=' + ','.join(sorted(protein_go[pid])))
            if pid in protein_ec:
                extra.append('ec_number=' + ','.join(sorted(protein_ec[pid])))
            if pid in protein_pfam:
                extra.append('Dbxref=' + ','.join('Pfam:' + p for p in sorted(protein_pfam[pid])))
            if pid in protein_product:
                extra.append('product=' + protein_product[pid].replace(';', ' '))
            if extra:
                fields[8] = fields[8].rstrip(';') + ';' + ';'.join(extra)
            fout.write('\\t'.join(fields) + '\\n')

    # Write summary
    n_proteins = sum(1 for line in open("${proteins}") if line.startswith('>'))
    n_annotated = len(annotations)
    n_total_ann = sum(len(v) for v in annotations.values())
    n_pfam = sum(1 for pid, anns in annotations.items() if any(a[0]=='Pfam' for a in anns))
    n_product = len(protein_product)
    with open("${sample_id}_summary.txt", 'w') as out:
        out.write(f"Sample: ${sample_id}\\n")
        out.write(f"Total proteins: {n_proteins}\\n")
        out.write(f"Annotated proteins: {n_annotated}\\n")
        out.write(f"Total annotations: {n_total_ann}\\n")
        out.write(f"With Pfam domains: {n_pfam}\\n")
        out.write(f"With product name: {n_product}\\n")
    """
}
