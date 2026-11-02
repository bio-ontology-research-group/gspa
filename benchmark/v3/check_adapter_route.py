"""Generic adapter route certification (row 3): genuine tool run, typed outcomes, route parity.

Driven by a per-adapter JSON spec (see benchmark/v3/results/adapter-certification-v1/*/spec.json).
For one adapter it runs the real GSPA CLI JAR through: annotate (genuine tool on PATH),
annotate+integration, direct integration of the raw claims, and then a set of injected cases
executed through a stub executable that replays variants of the genuine output
(header-only, malformed statistic, truncated row, dropped protein, unknown protein) plus a
missing-executable case. Every input is hashed before and after; the oracle is computed by
`adapter_oracles.py` without GSPA code. Adapter evidence only, not a biological benchmark.

Spec keys:
  adapter            predictor name as it appears in <genome>_predictor_status.json
  mode               "proteins" | "genome"
  input              {"proteins": faa} or {"genome": fna, "gff3": gff, "proteins": faa}
  config             GSPA config (JSON) for the genuine run
  path_dirs          directories holding the genuine tool(s), placed first on PATH
  tool               executable name the adapter resolves (stubbed in injected cases)
  stub_kind          how the stub finds its output location (see STUB)
  executable_config_key  JSON path to an executablePath field, or null (PATH-resolved tool)
  native_output      the genuine tool output file (used to derive injected variants)
  oracle_reader      name in adapter_oracles.READERS, or "regions" for genome mode
  genuine_outcome    expected outcome of the genuine run (SUCCESS | SUCCESS_ZERO | PARTIAL)
  cases              {case_name: {"expected": OUTCOME, "detail_contains": [..], "claims": "none"|"subset"|"same"}}
  bound              extra files/directories to hash
  genome_oracle      genome mode: {"gff3": path} for CDS attribution
  cli_threads        value passed to `annotate --threads` (default 1; it overrides the config's
                     per-predictor thread counts, so genome-wide search tools set it explicitly)
  timeout_seconds    wall-clock limit per CLI run (default 3600)
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adapter_oracles  # noqa: E402

STUB = r'''#!/usr/bin/python3
import os, shutil, sys
kind = os.environ['GSPA_STUB_KIND']; payload = os.environ.get('GSPA_STUB_PAYLOAD'); argv = sys.argv[1:]
def after(flag):
    return argv[argv.index(flag) + 1]
def place(dst, name='output'):
    if payload is None:
        return
    src = os.path.join(payload, name)
    if os.path.exists(src):
        os.makedirs(os.path.dirname(dst) or '.', exist_ok=True); shutil.copy(src, dst)
if argv[:1] in (['--version'], ['version'], ['-h'], ['--help'], ['-c']) or (kind == 'run_dbcan' and '--help' in argv):
    if kind == 'hmmsearch':  # the HMMER probe requires the real banner lines
        print('# hmmsearch :: search profile(s) against a sequence database\n# HMMER 3.4 (stub)')
    elif kind == 'run_dbcan':  # the dbCAN probe checks the flags of CAZyme_annotation --help
        print('stub --db_dir --input_raw_data --mode --threads')
    else:
        print('stub ' + kind + ' 0.0')
    sys.exit(0)
if kind == 'diamond':
    place(after('--out'))
elif kind == 'mmseqs':
    if argv[0] == 'convertalis':
        place(argv[4])
elif kind == 'foldseek':
    if argv[0] == 'easy-search':
        place(argv[3])
elif kind == 'minced':
    place(argv[-1])
elif kind == 'python-driver':
    place(argv[2]); place(os.path.join(os.path.dirname(argv[2]), 'coverage.tsv'), 'coverage')
elif kind == 'hmmsearch':
    place(after('--domtblout'))
elif kind == 'run_dbcan':
    place(os.path.join(after('--output_dir'), 'overview.tsv'))
elif kind == 'gapseq':
    place(os.path.join(after('-f'), 'query-Reactions.tbl'))
elif kind == 'genomic-sidecar':
    rows = open(after('--manifest')).read().splitlines()
    cols = rows[0].split('\t'); row = dict(zip(cols, rows[1].split('\t')))
    place(os.path.join(row['output_dir'], row['tag'] + '.' + after('--predictor') + '.genomic.tsv'))
else:
    sys.exit('unknown stub kind ' + kind)
sys.exit(0)
'''


CAPTURE = r'''#!/usr/bin/python3
"""Pass-through wrapper for run_genomic_predictors.py that keeps a copy of the sidecar TSV the
adapter is about to parse (the adapter deletes its scratch dir). Needed for tools that are not
deterministic between runs, so the oracle is read from the same run that produced the claims."""
import os, shutil, subprocess, sys
args = sys.argv[1:]
rows = open(args[args.index('--manifest') + 1]).read().splitlines()
row = dict(zip(rows[0].split('\t'), rows[1].split('\t')))
r = subprocess.run([sys.executable, os.environ['GSPA_CAPTURE_REAL']] + args)
src = os.path.join(row['output_dir'], row['tag'] + '.' + args[args.index('--predictor') + 1] + '.genomic.tsv')
if os.path.exists(src):
    shutil.copy(src, os.path.join(os.environ['GSPA_CAPTURE_DIR'], os.environ['GSPA_CAPTURE_NAME']))
sys.exit(r.returncode)
'''


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def tree(path):
    p = Path(path)
    return [p] if p.is_file() else sorted(q for q in p.rglob('*') if q.is_file())


def set_key(cfg, key, value):
    d = cfg
    for k in key[:-1]:
        d = d.setdefault(k, {})
    d[key[-1]] = value


KNOWN_CDS_IDS = set()


def cds_id(protein_id):
    """CDS-translated proteins are uniquified as <contig>:<protein_id>[.n] (gspa.io.CdsTranslator.uniqueId);
    the oracle keys by the GFF protein_id. A suffix is stripped only when the GFF knows the shorter id."""
    import re
    if protein_id in KNOWN_CDS_IDS:
        return protein_id
    base = re.sub(r'^[^:]+:', '', protein_id)
    if base in KNOWN_CDS_IDS:
        return base
    shorter = re.sub(r'\.\d+$', '', base)
    return shorter if shorter in KNOWN_CDS_IDS else base


def canonical(rows):
    return sorted(json.dumps(r, sort_keys=True) for r in rows)


def claim_view(row):
    """Reduce an exported raw claim to the oracle's comparison fields (identity, function, score, hit)."""
    m = row.get('metadata') or {}
    return {'protein_id': row['protein_id'], 'function_type': row['function_type'],
            'function_id': row['function_id'], 'raw_score': round(float(row['raw_score']), 6),
            'hit': m.get('hit', m.get('vfdbId'))}


oracle_view = claim_view


# ---- injected-output variants derived from the genuine native output -------------------------

def variants(spec, native, out):
    """Write one payload directory per injected case; return {case: payload_dir}."""
    text = native.read_text()
    lines = text.splitlines()
    kind = spec['stub_kind']
    payloads = {}

    def write(case, content, coverage=None):
        d = out / 'payload' / case
        d.mkdir(parents=True)
        (d / 'output').write_text(content)
        if coverage is not None:
            (d / 'coverage').write_text(coverage)
        payloads[case] = d

    header = [l for l in lines if l.startswith('#') or (kind == 'python-driver' and l.startswith('protein_id'))]
    data = [l for l in lines if l.strip() and l not in header]
    proteins = spec.get('input_protein_ids', [])
    full_cov = ''.join(f'{p}\tok\t\n' for p in proteins)
    full_cov = 'protein_id\tstatus\tmessage\n' + full_cov
    if kind in ('diamond', 'mmseqs', 'foldseek', 'minced'):
        write('header_only', '')  # documented: zero bytes when nothing aligned / no array
    elif kind == 'hmmsearch':
        write('header_only', '\n'.join(header) + '\n')
    elif kind == 'python-driver':
        write('header_only', 'protein_id\tregion_start\tregion_end\tmean_score\n', full_cov)
    elif kind in ('run_dbcan', 'genomic-sidecar', 'gapseq'):
        # comment lines (gapseq writes '# gapseq version ...') plus the named-column header row
        first = next(i for i, l in enumerate(lines) if l.strip() and not l.startswith('#'))
        header = lines[:first + 1]
        data = [l for l in lines[first + 1:] if l.strip()]
        write('header_only', '\n'.join(header) + '\n')
    if not data:
        return payloads  # genuine zero: only the header-only variant is derivable
    first = data[0]
    cols = first.split('\t') if kind != 'hmmsearch' else first.split()
    # malformed statistic: replace the numeric score column with text
    num_col = spec['malformed_column']
    bad = list(cols)
    bad[num_col] = 'not-a-number'
    joiner = '\t' if kind != 'hmmsearch' else ' '
    body = lambda rows: '\n'.join(header + rows) + '\n'
    write('malformed', body([joiner.join(bad)] + data[1:]), full_cov)
    # truncated row: drop the last two columns of the first data row (hmmsearch domtbl rows end
    # in a free-text description, so there the row is cut inside the fixed 22 columns)
    keep = cols[:15] if kind == 'hmmsearch' else cols[:-2]
    write('truncated_row', body([joiner.join(keep)] + data[1:]), full_cov)
    # dropped protein: remove every row of the protein named in the spec
    drop = spec['dropped_protein']
    kept = [l for l in data if not l.split()[0] == drop]
    dropped_cov = 'protein_id\tstatus\tmessage\n' + ''.join(f'{p}\tok\t\n' for p in proteins if p != drop)
    write('dropped_protein', body(kept), dropped_cov)
    # unknown protein: add a row naming a protein that was never submitted
    ghost = list(cols)
    ghost[0] = 'ghost_protein_not_in_input'
    write('unknown_protein', body(data + [joiner.join(ghost)]), full_cov + 'ghost_protein_not_in_input\tok\t\n')
    return payloads


# ---- main -----------------------------------------------------------------------------------

def run(spec_path, jar, java, out):
    spec = json.loads(Path(spec_path).read_text())
    out.mkdir(parents=True, exist_ok=False)
    base = Path(spec_path).resolve().parent
    rel = lambda p: (base / p) if not Path(p).is_absolute() else Path(p)
    native = rel(spec['native_output'])
    inputs = {k: rel(v) for k, v in spec['input'].items()}
    # protein ids for coverage payloads
    fasta = inputs.get('proteins')
    ids = [l[1:].split()[0] for l in fasta.read_text().splitlines() if l.startswith('>')] if fasta else []
    spec['input_protein_ids'] = ids

    bound = {}
    for p in [jar, java, Path(__file__), Path(adapter_oracles.__file__), Path(spec_path), native, *inputs.values()]:
        bound[str(p)] = sha(p)
    for extra in spec.get('bound', []):
        for p in tree(rel(extra)):
            bound[str(p)] = sha(p)
    for d in spec.get('path_dirs', []):
        for p in tree(d):
            if p.stat().st_size < (1 << 30) and os.access(p, os.R_OK):
                bound[str(p)] = sha(p)

    # oracle
    cds = None
    if spec['oracle_reader'] == 'none':
        oracle = {'note': 'no claim oracle: the route under test must produce zero claims'}
    elif spec['oracle_reader'] in ('regions', 'sidecar_regions'):
        regions = (adapter_oracles.minced_regions if spec['oracle_reader'] == 'regions' else adapter_oracles.sidecar_regions)(native)
        cds = adapter_oracles.gff_cds(rel(spec['genome_oracle']['gff3']))
        KNOWN_CDS_IDS.update(c['protein_id'] for c in cds)
        expected_proteins = adapter_oracles.region_attribution(regions, cds)
        oracle = {'regions': regions, 'attributed_proteins': expected_proteins}
    else:
        oracle = adapter_oracles.READERS[spec['oracle_reader']](native)
    (out / 'oracle.json').write_text(json.dumps(oracle, indent=1) + '\n')

    # stub executable
    stubbin = out / 'stubbin'
    stubbin.mkdir()
    stub = stubbin / spec['tool']
    stub.write_text(STUB)
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    emptybin = out / 'emptybin'
    emptybin.mkdir()
    for helper in ('bash', 'sh', 'env'):
        os.symlink('/usr/bin/' + helper, emptybin / helper)
        os.symlink('/usr/bin/' + helper, stubbin / helper)
    payloads = variants(spec, native, out)

    workflow = out / 'workflow.json'
    workflow.write_text(json.dumps({'lite': True, 'enabledPriors': 'none', 'genomeId': 'input'}) + '\n')
    bound[str(workflow)] = sha(workflow)
    tmp = out / 'tmp'
    tmp.mkdir()
    runs, statuses, claims, checks, observations = [], {}, {}, [], []

    def check(name, ok, detail=''):
        checks.append({'check': name, 'ok': bool(ok), 'detail': detail})
        if not ok:
            print('FAIL', name, detail)

    def execute(name, command, path_dirs, stub_env=None):
        env = dict(os.environ, TMPDIR=str(tmp), PATH=':'.join(str(p) for p in path_dirs))
        if stub_env:
            env.update(stub_env)
        with (out / (name + '.log')).open('w') as log:
            r = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                               timeout=spec.get('timeout_seconds', 3600))
        runs.append({'name': name, 'command': command, 'exit_code': r.returncode, 'path': env['PATH'],
                     'stub': {k: v for k, v in (stub_env or {}).items()}})
        (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        return r.returncode

    def annotate(name, cfg, path_dirs, integrate=False, stub_env=None):
        cfg_file = out / (name + '.config.json')
        cfg_file.write_text(json.dumps(cfg, indent=2) + '\n')
        folder = out / name
        cmd = [str(java), '-XX:ActiveProcessorCount=4', f'-Djava.io.tmpdir={tmp}', '-jar', str(jar), 'annotate']
        if spec['mode'] == 'genome':
            cmd += ['--genome', str(inputs['genome']), '--gff3', str(inputs['gff3'])]
            if 'proteins' in inputs:  # otherwise the CLI translates the GFF CDS itself
                cmd += ['--proteins', str(inputs['proteins'])]
            stem = inputs['genome'].name.rsplit('.', 1)[0]
        else:
            cmd += ['--proteins', str(inputs['proteins'])]
            stem = inputs['proteins'].name.rsplit('.', 1)[0]
        # --threads overrides the per-predictor thread counts of the config (AnnotateCommand), so a
        # genome-wide search tool (gapseq find -M nucl: 155,682 reference sequences by tblastn)
        # needs the spec's cli_threads; 1 is enough for the protein-fixture adapters.
        cmd += ['--config', str(cfg_file), '--output', str(folder), '--claims-out', str(folder / 'raw.jsonl'),
                '--threads', str(spec.get('cli_threads', 1)), '--no-metrics']
        if integrate:
            cmd += ['--integration-config', str(workflow)]
        rc = execute(name, cmd, path_dirs, stub_env)
        check(name + ':exit0', rc == 0, f'exit {rc}')
        status_file = folder / f'{stem}_predictor_status.json'
        st = json.loads(status_file.read_text())['predictors'] if status_file.is_file() else []
        mine = [s for s in st if s['predictor'] == spec['adapter']]
        statuses[name] = mine[0] if mine else {'predictor': spec['adapter'], 'outcome': 'ABSENT', 'detail': 'not in manifest'}
        raw = folder / 'raw.jsonl'
        rows = [json.loads(l) for l in raw.read_text().splitlines()] if raw.is_file() else []
        claims[name] = [r for r in rows if r.get('source') == spec['adapter']]
        return folder, stem

    cfg = spec['config']
    genuine_dirs = [Path(p) for p in spec['path_dirs']] + [Path('/usr/bin')]
    capture_env = None
    if spec.get('capture_native'):   # nondeterministic genomic tools: oracle from the same run
        wrapper = out / 'capture_sidecar.py'
        wrapper.write_text(CAPTURE)
        (out / 'captured').mkdir()
        real = cfg
        for k in spec['executable_config_key']:
            real = real[k]
        cfg = json.loads(json.dumps(cfg))
        set_key(cfg, spec['executable_config_key'], str(wrapper))
        capture_env = {'GSPA_CAPTURE_REAL': real, 'GSPA_CAPTURE_DIR': str(out / 'captured')}

    # 1. genuine annotation, integrated annotation, direct integration
    folder, stem = annotate('annotate', cfg, genuine_dirs,
                            stub_env=dict(capture_env, GSPA_CAPTURE_NAME='annotate.native.tsv') if capture_env else None)
    if capture_env and (out / 'captured' / 'annotate.native.tsv').is_file():
        live = adapter_oracles.sidecar_regions(out / 'captured' / 'annotate.native.tsv')
        observations.append({'observation': 'direct_run_vs_annotate_run_regions',
                             'identical': sorted((r['contig'], r['start'], r['end']) for r in live) == sorted((r['contig'], r['start'], r['end']) for r in oracle['regions']),
                             'direct_run': [(r['contig'], r['start'], r['end']) for r in oracle['regions']],
                             'annotate_run': [(r['contig'], r['start'], r['end']) for r in live]})
        oracle = {'regions': live, 'attributed_proteins': adapter_oracles.region_attribution(live, cds)}
        (out / 'oracle.annotate-run.json').write_text(json.dumps(oracle, indent=1) + '\n')
    check('annotate:outcome', statuses['annotate']['outcome'] == spec['genuine_outcome'],
          f"{statuses['annotate']['outcome']} != {spec['genuine_outcome']} ({statuses['annotate'].get('detail')})")
    for needle in spec.get('genuine_detail_contains', []):
        check(f'annotate:detail~{needle}', needle in (statuses['annotate'].get('detail') or ''), statuses['annotate'].get('detail'))
    if spec.get('genuine_claims') == 'none':
        check('annotate:no_claims', not claims['annotate'], f"{len(claims['annotate'])} claims")
    if spec['oracle_reader'] == 'none':
        pass
    elif spec['oracle_reader'] in ('regions', 'sidecar_regions'):
        got = sorted({cds_id(r['protein_id']) for r in claims['annotate']})
        known = spec.get('attribution_known_gaps', {})   # {protein_id: reason}: documented, not excused silently
        missing = sorted(set(oracle['attributed_proteins']) - set(got))
        extra = sorted(set(got) - set(oracle['attributed_proteins']))
        check('annotate:attribution', not extra and all(m in known for m in missing),
              f'{len(got)} attributed vs {len(oracle["attributed_proteins"])} expected; missing={missing} extra={extra}')
        for m in missing:
            observations.append({'observation': 'attribution_gap', 'protein_id': m, 'reason': known.get(m)})
        gff3 = folder / f'{stem}.gff3'
        region_rows = [l for l in gff3.read_text().splitlines() if not l.startswith('#') and l.count('\t') >= 8
                       and l.split('\t')[1] == spec['adapter'] and l.split('\t')[2] in ('region', 'repeat_region')] if gff3.is_file() else []
        got_regions = sorted((l.split('\t')[0], int(l.split('\t')[3]), int(l.split('\t')[4])) for l in region_rows)
        exp_regions = sorted((r['contig'], r['start'], r['end']) for r in oracle['regions'])
        check('annotate:regions_retained', got_regions == exp_regions, f'{got_regions} vs {exp_regions}')
    else:
        got = canonical(claim_view(r) for r in claims['annotate'])
        exp = canonical(oracle_view(r) for r in oracle)
        check('annotate:claims=oracle', got == exp, f'{len(got)} exported vs {len(exp)} expected; '
              f'missing={[e for e in exp if e not in got][:3]} extra={[g for g in got if g not in exp][:3]}')
        for row in claims['annotate']:
            m = row.get('metadata') or {}
            exp_meta = [o['metadata'] for o in oracle if oracle_view(o) == claim_view(row)]
            if exp_meta:
                sub = {k: v for k, v in exp_meta[0].items() if k in m}
                check(f'annotate:metadata:{row["protein_id"]}:{row["function_id"]}',
                      all(abs(float(m[k]) - float(v)) < 1e-6 if isinstance(v, (int, float)) else str(m[k]) == str(v)
                          for k, v in sub.items()), f'{sub} vs {m}')
    for needle in spec.get('genuine_skipped_contains', []):
        check(f'annotate:skipped~{needle}', needle in (statuses['annotate'].get('skipped_inputs') or []), str(statuses['annotate'].get('skipped_inputs')))
    annotate('annotate_integrated', cfg, genuine_dirs, integrate=True,
             stub_env=dict(capture_env, GSPA_CAPTURE_NAME='annotate_integrated.native.tsv') if capture_env else None)
    same = canonical(claims['annotate']) == canonical(claims['annotate_integrated'])
    if spec.get('nondeterministic'):
        # the tool itself may differ between two genuine runs; record it, and check the route
        # parity on the integrated run's own raw claims
        observations.append({'observation': 'raw_claims_between_two_genuine_runs', 'identical': same,
                             'annotate': len(claims['annotate']), 'annotate_integrated': len(claims['annotate_integrated'])})
        parity_source = out / 'annotate_integrated' / 'raw.jsonl'
    else:
        check('parity:raw_claims', same)
        parity_source = out / 'annotate' / 'raw.jsonl'
    direct = out / 'direct.tsv'
    rc = execute('direct_integrate', [str(java), f'-Djava.io.tmpdir={tmp}', '-jar', str(jar), 'integrate', '--claims',
                 str(parity_source), '--workflow-config', str(workflow), '--out', str(direct)], genuine_dirs)
    check('direct_integrate:exit0', rc == 0)
    integrated = out / 'annotate_integrated' / f'{stem}_integrated.tsv'
    check('parity:integration_bytes', direct.is_file() and integrated.is_file() and direct.read_bytes() == integrated.read_bytes())

    # 2. injected cases through the stub executable
    for case, exp in spec['cases'].items():
        c = json.loads(json.dumps(spec['config']))
        if case == 'unavailable':
            if spec.get('executable_config_key'):
                set_key(c, spec['executable_config_key'], str(out / 'absent-tool'))
            annotate(case, c, [emptybin])
        else:
            if 'output' in exp:  # spec-supplied payload (e.g. attribution cases for genome tools)
                d = out / 'payload' / case
                d.mkdir(parents=True)
                (d / 'output').write_text(exp['output'])
                if 'coverage' in exp:
                    (d / 'coverage').write_text(exp['coverage'])
                payloads[case] = d
            if case not in payloads:
                check(case + ':payload', False, 'no payload derivable (genuine output had no data rows)')
                continue
            if spec.get('executable_config_key'):
                set_key(c, spec['executable_config_key'], str(stub))
            annotate(case, c, [stubbin], stub_env={'GSPA_STUB_KIND': spec['stub_kind'], 'GSPA_STUB_PAYLOAD': str(payloads[case])})
        st = statuses[case]
        check(f'{case}:outcome', st['outcome'] == exp['expected'], f"{st['outcome']} != {exp['expected']} ({st.get('detail')})")
        for needle in exp.get('detail_contains', []):
            check(f'{case}:detail~{needle}', needle in (st.get('detail') or ''), st.get('detail'))
        if exp.get('claims') == 'none':
            check(f'{case}:no_claims', not claims[case], f'{len(claims[case])} claims')
        elif exp.get('claims') == 'without_dropped':
            dropped = spec['dropped_protein']
            check(f'{case}:dropped_absent', all(r['protein_id'] != dropped for r in claims[case]))
            check(f'{case}:others_kept', canonical(claim_view(r) for r in claims[case]) ==
                  canonical(claim_view(r) for r in claims['annotate'] if r['protein_id'] != dropped))
        if exp.get('skipped_contains'):
            check(f'{case}:skipped_inputs', exp['skipped_contains'] in (st.get('skipped_inputs') or []), str(st.get('skipped_inputs')))
        if 'expected_proteins' in exp:
            got = sorted({cds_id(r['protein_id']) for r in claims[case]})
            check(f'{case}:attributed_proteins', got == sorted(exp['expected_proteins']), f'{got} vs {exp["expected_proteins"]}')

    # 3. inputs unchanged
    for p, digest in bound.items():
        check('unchanged:' + p, sha(p) == digest)
    ok = all(c['ok'] for c in checks)
    report = {'schema': 'gspa.adapter-route-replay/1', 'adapter': spec['adapter'],
              'status': 'verified' if ok else 'failed', 'checks_total': len(checks),
              'checks_failed': [c for c in checks if not c['ok']], 'observations': observations, 'predictor_outcomes': statuses,
              'claims_exported': {k: len(v) for k, v in claims.items()}, 'inputs': bound,
              'outputs': {str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file() and 'tmp' not in p.parts}}
    (out / 'checks.json').write_text(json.dumps(checks, indent=1) + '\n')
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'adapter': spec['adapter'], 'status': report['status'], 'checks': len(checks),
                      'failed': len(report['checks_failed']), 'outcomes': {k: v['outcome'] for k, v in statuses.items()}}))
    return ok


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--spec', type=Path, required=True)
    ap.add_argument('--jar', type=Path, required=True)
    ap.add_argument('--java', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    sys.exit(0 if run(a.spec.resolve(), a.jar.resolve(), a.java.resolve(), a.out.resolve()) else 1)
