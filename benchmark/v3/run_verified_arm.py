"""Run one planned WP4 arm through ``verified_batch.run_batch`` (receipt ``verified_run.json``).

The WP4 planner (``prepare_wp4_batches.py plan``) writes each run as ``command.json`` (an argv) plus a batch
manifest whose ``out_path`` rows lie directly inside the run's ``inference/`` directory. This wrapper executes that
argv under the verified-batch contract: existing outputs are archived to ``history/``, the command is rerun, and
the receipt binds the command, JAR, every input file, the runner sources and every prediction and stage ledger.
A learned run may use the documented unsupported-learned fallback (as ``check_source_arm_context_models.py``),
which ``run_batch`` validates per genome; a run without ``--learned-model`` never may.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import verified_batch as V


def load_command(path):
    command = json.loads(Path(path).read_text())
    if not (isinstance(command, list) and command and all(isinstance(t, str) for t in command)):
        raise ValueError(f"command.json is not a string argv: {path}")
    return command


def run(command_json, manifest, out_dir, timeout=7200):
    command = load_command(command_json)
    learned = "--learned-model" in command
    sources = [Path(__file__), Path(command_json)]
    return V.run_batch(command, Path(manifest), Path(out_dir), sources=sources, timeout=timeout,
                       allow_unsupported_learned=learned)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--command-json", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=7200)
    args = parser.parse_args(argv)
    receipt = run(args.command_json, args.manifest, args.out_dir, args.timeout)
    print(json.dumps({"returncode": receipt["returncode"], "artifacts": len(receipt["artifacts"]),
                      "required_stages": receipt["spec"]["required_stages"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
