"""Independent readers of genuine tool outputs → expected GSPA raw claims.

Each reader re-implements the documented output format of the external tool and the
adapter's documented scoring rule without importing or calling GSPA code. The result is a
list of claim dicts (protein_id, function_type, function_id, raw_score, metadata subset)
that `check_adapter_route.py` compares with the raw claims the CLI actually exports.
This is adapter evidence (format, dispatch, attribution), not a biological benchmark.
"""
import json
import re
from pathlib import Path

GO = re.compile(r'GO:\d{7}')
EC = re.compile(r'EC:[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+')


def _rows(path, min_fields):
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < min_fields:
            raise ValueError('short row: ' + line[:80])
        yield fields


def _go_ec(pid, header, score, hit, extra):
    out = []
    for go in GO.findall(header):
        out.append({'protein_id': pid, 'function_type': 'GO', 'function_id': go, 'raw_score': score,
                    'metadata': dict(hit=hit, **extra)})
    for ec in EC.findall(header):
        out.append({'protein_id': pid, 'function_type': 'EC', 'function_id': ec, 'raw_score': score,
                    'metadata': {'hit': hit}})
    return out


def diamond(path):
    """qseqid sseqid pident length qlen slen evalue bitscore stitle; score = min(1, pident/100)."""
    out = []
    for f in _rows(path, 9):
        pident = float(f[2])
        out += _go_ec(f[0], f[8], min(1.0, pident / 100.0), f[1], {'pident': pident, 'evalue': float(f[6])})
    return out


def mmseqs2(path):
    """query,target,pident,alnlen,qlen,tlen,evalue,bits,theader; pident already in percent."""
    out = []
    for f in _rows(path, 9):
        pident = float(f[2])
        out += _go_ec(f[0], f[8], min(1.0, pident / 100.0), f[1], {'pident': pident, 'evalue': float(f[6])})
    return out


def foldseek_prostt5(path):
    """query,target,pident,evalue,bits,theader (ProstT5 mode: no TM-score, score = pident/100)."""
    out = []
    for f in _rows(path, 6):
        pident = float(f[2])
        out += _go_ec(f[0], f[5], min(1.0, pident / 100.0), f[1], {'pident': pident, 'evalue': float(f[3])})
    return out


def vfdb(path):
    """qseqid sseqid pident length evalue bitscore stitle; value = text after the first '(gene)' token."""
    out = []
    for f in _rows(path, 7):
        pident = float(f[2])
        m = re.search(r'\((\w+)\)\s+(.*)', f[6])
        value = (m.group(2) if m else f[6][:100]).strip()
        out.append({'protein_id': f[0], 'function_type': 'VFDB', 'function_id': value,
                    'raw_score': min(1.0, pident / 100.0),
                    'metadata': {'vfdbId': f[1], 'pident': pident, 'evalue': float(f[4])}})
    return out


def metapredict(path, min_len=10, min_score=0.5):
    """protein_id region_start region_end mean_score (1-based inclusive); adapter filters length/score."""
    out = []
    for line in Path(path).read_text().splitlines()[1:]:
        if not line.strip():
            continue
        pid, s, e, score = line.split('\t')[:4]
        s, e, score = int(s), int(e), float(score)
        if e - s + 1 < min_len or score < min_score:
            continue
        out.append({'protein_id': pid, 'function_type': 'DISORDER', 'function_id': 'disorder',
                    'raw_score': score, 'metadata': {}})
    return out


def antifam_domtbl(path):
    """hmmsearch --domtblout; one PSEUDOGENE quality observation per (protein, model), score 1.0 under --cut_ga."""
    seen = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        acc = re.sub(r'\.\d+$', '', f[4])
        key = (f[0], acc)
        if key in seen:
            seen[key]['metadata']['domainCount'] += 1
            continue
        seen[key] = {'protein_id': f[0], 'function_type': 'PSEUDOGENE', 'function_id': acc, 'raw_score': 1.0,
                     'metadata': {'name': f[3], 'modelAccession': f[4], 'domainCount': 1,
                                  'seqEvalue': float(f[6]), 'domainEvalue': float(f[12])}}
    return list(seen.values())


def minced_regions(path):
    """MinCED -gffFull: repeat_region rows are arrays (contig, start, end, 1-based inclusive)."""
    regions = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f = line.split('\t')
        if len(f) >= 9 and f[2] == 'repeat_region':
            regions.append({'contig': f[0], 'start': int(f[3]), 'end': int(f[4])})
    return regions


def sidecar_regions(path):
    """GSPA genomic sidecar TSV (contig_id, region_start, region_end, region_type, score, attributes)."""
    regions = []
    for line in Path(path).read_text().splitlines()[1:]:
        if not line.strip():
            continue
        f = line.split('\t')
        regions.append({'contig': f[0], 'start': int(f[1]), 'end': int(f[2]), 'type': f[3]})
    return regions


def phispy_coordinates(path):
    """PhiSpy prophage_coordinates.tsv: pp_number, contig, start, end, ..."""
    return [{'contig': f[1], 'start': int(f[2]), 'end': int(f[3])} for f in
            (l.split('\t') for l in Path(path).read_text().splitlines() if l.strip())]


def gff_cds(path):
    """RefSeq GFF3 CDS rows → (contig, start, end, protein_id) using the Name/protein_id attribute."""
    cds = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        attrs = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        pid = attrs.get('protein_id') or attrs.get('Name')
        if pid:
            cds.append({'contig': f[0], 'start': int(f[3]), 'end': int(f[4]), 'protein_id': pid})
    return cds


def region_attribution(regions, cds):
    """Protein ids whose CDS overlaps a region (same contig, 1-based inclusive overlap)."""
    hits = set()
    for r in regions:
        for c in cds:
            if c['contig'] == r['contig'] and c['start'] <= r['end'] and r['start'] <= c['end']:
                hits.add(c['protein_id'])
    return sorted(hits)


FAMILY = re.compile(r'^((?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?)(?:_e\d+)?(?:\((\d+)-(\d+)\))?$')


def dbcan_overview(path):
    """run_dbcan 5 overview.tsv: one unscored CAZYME claim per (protein, CAZy family) over the
    three method columns; dbCAN_sub subfamilies (_eN) collapse to their parent family."""
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    header = lines[0].split('\t')
    cols = {h: i for i, h in enumerate(header)}
    out = []
    for line in lines[1:]:
        f = line.split('\t')
        families = []
        for col in ('dbCAN_hmm', 'HMMER', 'dbCAN', 'dbCAN_sub', 'eCAMI', 'DIAMOND'):
            if col not in cols or f[cols[col]] == '-':
                continue
            for token in f[cols[col]].split('+'):
                m = FAMILY.match(token.strip())
                if not m:
                    raise ValueError('bad family token ' + token)
                if m.group(1) not in families:
                    families.append(m.group(1))
        for fam in families:
            out.append({'protein_id': f[cols['Gene ID']], 'function_type': 'CAZYME', 'function_id': fam,
                        'raw_score': 0.0, 'metadata': {'proteinToolCount': int(f[cols['#ofTools']])}})
    return out


READERS = {'dbcan': dbcan_overview, 'diamond': diamond, 'mmseqs2': mmseqs2, 'foldseek': foldseek_prostt5, 'vfdb': vfdb,
           'metapredict': metapredict, 'antifam': antifam_domtbl}


if __name__ == '__main__':
    import sys
    print(json.dumps(READERS[sys.argv[1]](sys.argv[2]), indent=1))
