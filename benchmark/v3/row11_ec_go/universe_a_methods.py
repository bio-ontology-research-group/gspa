#!/usr/bin/env python3
"""Universe A deterministic/learned methods and agent driver using EngineConsistentGuardrail."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_stage3_agent as E  # noqa: E402
import v3_runner as runner  # noqa: E402
from row11_ec_go import guardrail as G  # noqa: E402
from row11_ec_go.guardrail import EngineConsistentGuardrail, load as load_guardrail  # noqa: E402

MODES = E.MODES


def suggest(public: dict, mode: str, guardrail: EngineConsistentGuardrail, solve=None) -> dict:
    return E.suggest(public, mode, guardrail, solve=solve)


def suggest_agent(
    public: dict,
    guardrail: EngineConsistentGuardrail,
    model: str = "qwen3.8-27b",
    api_base: str | None = None,
    max_turns: int = 2,
    workers: int = 4,
) -> dict:
    return E.suggest_agent(public, guardrail, model, api_base, max_turns, workers)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="phase", required=True)
    for phase in ("suggest", "agent"):
        cmd = sub.add_parser(phase)
        cmd.add_argument("--public", required=True, type=Path)
        cmd.add_argument("--go-obo", required=True, type=Path)
        if phase == "suggest":
            cmd.add_argument("--mode", required=True, choices=MODES)
        else:
            cmd.add_argument("--model", default="qwen3.8-27b")
            cmd.add_argument("--api-base", default=None)
            cmd.add_argument("--max-turns", type=int, default=2)
            cmd.add_argument("--workers", type=int, default=4)
        cmd.add_argument("--out", required=True, type=Path)

    args = parser.parse_args(argv)
    code_paths = (
        Path(__file__).resolve(),
        Path(G.__file__).resolve(),
        Path(E.__file__).resolve(),
        Path(E.agent.__file__).resolve(),
    )
    code_identity = {str(path.relative_to(ROOT)): runner.sha256_file(str(path)) for path in code_paths}
    resources = ROOT / "gspa-core" / "src" / "main" / "resources"
    guard = load_guardrail(resources, args.go_obo)

    if args.phase == "suggest":
        result = E.suggest(E.read(args.public), args.mode, guard)
    else:
        result = E.suggest_agent(
            E.read(args.public),
            guard,
            args.model,
            args.api_base,
            args.max_turns,
            args.workers,
        )

    result["decision_code_sha256"] = code_identity
    result["constraint_inputs"] = [runner.hash_record(str(args.go_obo))] + [
        runner.hash_record(str(path))
        for directory in ("taxon-constraints", "anatomical-constraints")
        for path in sorted((resources / directory).glob("*.tsv"))
    ]
    if code_identity != {str(path.relative_to(ROOT)): runner.sha256_file(str(path)) for path in code_paths}:
        raise RuntimeError("decision code changed during execution")
    E.write(args.out, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
