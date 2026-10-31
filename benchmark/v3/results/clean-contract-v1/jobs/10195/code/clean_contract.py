"""Typed CLEAN max-separation adapter. Native distances are never probabilities."""
from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import uuid

HEADER = ['protein_id', 'term', 'score', 'annotation_type', 'metadata_json']
KINDS = ('gmm_confidence', 'distance')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_mode(driver, kind):
    if kind not in KINDS:
        raise ValueError('CLEAN requires explicit --clean-score-kind gmm_confidence or distance')
    tree = ast.parse(driver.read_text())
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Name) and n.func.id == 'infer_maxsep']
    if len(calls) != 1 or len(calls[0].args) > 2 or any(k.arg is None for k in calls[0].keywords):
        raise ValueError('unsupported CLEAN driver: expected one explicit infer_maxsep call')
    argument = next((k.value for k in calls[0].keywords if k.arg == 'gmm'), ast.Constant(None))
    if not isinstance(argument, ast.Constant) or (argument.value is not None and not isinstance(argument.value, str)):
        raise ValueError('unsupported dynamic CLEAN gmm argument')
    observed = 'distance' if argument.value is None else 'gmm_confidence'
    if observed != kind:
        raise ValueError(f'declared CLEAN score kind {kind} conflicts with driver {observed}')
    return {'driver_sha256': sha(driver), 'driver_score_kind': observed,
            'gmm_path': argument.value}


def read_proteins(path):
    records, seen, current, sequence = [], set(), None, []
    def finish():
        if current is not None:
            value = ''.join(sequence)
            if not value:raise ValueError('empty protein sequence: ' + current)
            records.append((current, value))
    for line in Path(path).read_text().splitlines():
        if not line.strip():continue
        if line.startswith('>'):
            finish()
            tokens = line[1:].split()
            if not tokens or tokens[0] in seen:raise ValueError('empty or duplicated FASTA protein ID')
            current = tokens[0];seen.add(current);sequence = []
        else:
            if current is None:raise ValueError('sequence before FASTA header')
            sequence.append(''.join(line.split()))
    finish()
    if not records:raise ValueError('CLEAN requires a nonempty protein FASTA')
    return records


def rewrite(src, dst, min_score, *, score_kind, protein_ids, provenance):
    """Validate every native cell before publishing TSV+receipt; preserve all IDs."""
    src, dst = Path(src), Path(dst)
    receipt_path = dst.with_suffix('.receipt.json')
    dst.unlink(missing_ok=True);receipt_path.unlink(missing_ok=True)
    if score_kind not in KINDS:raise ValueError('CLEAN score kind is required')
    if not math.isfinite(min_score) or not 0 <= min_score <= 1:raise ValueError('invalid CLEAN minimum confidence')
    if score_kind == 'distance' and min_score != 0:
        raise ValueError('min-score is not a distance cutoff; use zero to retain all selected distances')
    if not protein_ids or len(set(protein_ids.values())) != len(protein_ids):
        raise ValueError('requires a complete one-to-one native/input protein map')
    selected, emitted, seen = 0, 0, set()
    common = {**provenance, 'native_output_sha256': sha(src), 'score_kind': score_kind,
              'score_calibrated': False, 'score_available': score_kind == 'gmm_confidence',
              'predictor_contract': 'clean:maxsep/1:' + score_kind + ':' + provenance['driver_sha256'],
              'candidate_coverage': 'native_maxsep_selected_only; not the complete EC distance matrix'}
    dst.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.clean-', dir=dst.parent) as work:
        temporary = Path(work) / 'output.tsv'
        with src.open(newline='') as source, temporary.open('w') as output:
            output.write('\t'.join(HEADER) + '\n')
            for line, fields in enumerate(csv.reader(source), 1):
                if len(fields) < 2 or fields[0] not in protein_ids or fields[0] in seen:
                    raise ValueError(f'CLEAN row {line}: missing predictions, unknown or duplicated protein ID')
                native_id = fields[0];seen.add(native_id);functions = set()
                for token in fields[1:]:
                    parts = token.strip().split('/')
                    if len(parts) != 2 or not re.fullmatch(r'EC:[1-7](?:\.(?:\d+|-)){3}', parts[0]):
                        raise ValueError(f'CLEAN row {line}: invalid EC/distance-or-confidence token')
                    label, raw = parts
                    if label in functions:raise ValueError(f'CLEAN row {line}: duplicated EC call')
                    functions.add(label)
                    native_score = float(raw)
                    if not math.isfinite(native_score) or native_score < 0 or (score_kind == 'gmm_confidence' and native_score > 1):
                        raise ValueError(f'CLEAN row {line}: invalid {score_kind}')
                    selected += 1
                    if score_kind == 'gmm_confidence' and native_score < min_score:continue
                    # A primitive placeholder is deliberately marked unavailable for distance calls.
                    # Native distances/text remain in metadata, with their lower-is-better direction.
                    score = native_score if score_kind == 'gmm_confidence' else 0.0
                    metadata = {**common, 'native_score': native_score, 'native_score_text': raw,
                                'native_protein_id': native_id,
                                'native_score_direction': 'higher' if score_kind == 'gmm_confidence' else 'lower'}
                    output.write('\t'.join([protein_ids[native_id], label, repr(score), 'EC',
                                            json.dumps(metadata, sort_keys=True, allow_nan=False)]) + '\n')
                    emitted += 1
        if seen != set(protein_ids):raise ValueError('CLEAN native output is incomplete for the input proteins')
        receipt = {**common, 'schema': 'gspa.clean-output/1', 'input_protein_ids': list(protein_ids.values()),
                   'native_to_input_ids': dict(protein_ids),
                   'input_proteins': len(protein_ids), 'native_proteins': len(seen),
                   'native_selected_calls': selected, 'emitted_calls': emitted,
                   'filtered_calls': selected - emitted, 'min_score': min_score,
                   'output_sha256': sha(temporary)}
        staged_receipt = Path(work) / 'receipt.json'
        staged_receipt.write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
        try:
            temporary.replace(dst);staged_receipt.replace(receipt_path)
        except Exception:
            dst.unlink(missing_ok=True);receipt_path.unlink(missing_ok=True);raise
    return receipt


def run(rows, args):
    if not args.model_dir:raise ValueError('CLEAN requires --model-dir')
    app = Path(args.model_dir).resolve()
    driver = app / 'CLEAN_infer_fasta.py'
    mode = validate_mode(driver, getattr(args, 'clean_score_kind', None))
    for row in rows:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', row.tag) or row.tag in ('.', '..'):
            raise ValueError('invalid CLEAN output tag')
        row.output_dir.mkdir(parents=True, exist_ok=True)
        output = row.output_dir / (row.tag + '.clean.tsv')
        output.unlink(missing_ok=True);output.with_suffix('.receipt.json').unlink(missing_ok=True)
        input_sha = sha(row.fasta_path)
        proteins = read_proteins(row.fasta_path)
        namespace = 'gspa_' + uuid.uuid4().hex
        identifiers = {f'{namespace}_p{i:08d}': pid for i, (pid, _) in enumerate(proteins)}
        stage = app / 'data/inputs' / (namespace + '.fasta')
        produced = app / 'results/inputs' / (namespace + '_maxsep.csv')
        stage.parent.mkdir(parents=True, exist_ok=True);produced.parent.mkdir(parents=True, exist_ok=True)
        owned = [stage, stage.with_suffix('.csv'), produced]
        owned.extend(app / 'data/esm_data' / (pid + '.pt') for pid in identifiers)
        try:
            stage.write_text(''.join('>' + native_id + '\n' + seq + '\n'
                                     for native_id, (_, seq) in zip(identifiers, proteins)))
            subprocess.run([sys.executable, str(driver), '--fasta_data', namespace], check=True, cwd=app)
            if not produced.is_file():raise ValueError('CLEAN did not produce this run\'s native output')
            if sha(row.fasta_path) != input_sha:raise ValueError('input FASTA changed during CLEAN execution')
            native_copy = row.output_dir / (row.tag + '.clean.native.csv')
            shutil.copyfile(produced, native_copy)
            rewrite(native_copy, output, args.min_score, score_kind=args.clean_score_kind,
                    protein_ids=identifiers, provenance={**mode, 'input_fasta_sha256': input_sha})
        finally:
            for path in owned:path.unlink(missing_ok=True)
