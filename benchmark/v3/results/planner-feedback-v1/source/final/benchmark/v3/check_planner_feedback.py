#!/usr/bin/env python3
"""Paired planner plumbing diagnostic. Synthetic network, no accuracy claim.

Run on a compute node. --live-qwen uses the configured unimatrix endpoint and
retains prompts, responses and usage, never authentication headers.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "scripts/tests"), str(Path(__file__).parent)]
import metabolic_experiments as ME
from test_metabolic_experiments import manifest, write_sbml
from run_planner_benchmark import UsageChat


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--before", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--live-qwen", action="store_true")
    p.add_argument("--key-file", type=Path)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    if args.key_file:
        os.environ["BORG_LLM_API_KEY"] = args.key_file.read_text().strip()
    spec = importlib.util.spec_from_file_location("before_feedback", args.before)
    before = importlib.util.module_from_spec(spec); sys.modules[spec.name] = before
    spec.loader.exec_module(before)
    model_path = args.out / "toy.xml"
    write_sbml(model_path)
    model = ME.load_sbml_model(model_path)
    old_model = before.load_sbml_model(model_path)

    def save(name, data):
        (args.out / name).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")

    summary = {"scope": "synthetic mechanism and routing; no biological or agent-benefit claim",
               "non_agent_pairs": [], "agent_runs": []}
    for strategy in ("fixed", "random", "evidence"):
        for seed in (0, 7):
            m = manifest(strategy=strategy, seed=seed, max_experiments=20, max_solver_calls=100)
            m["model"] = str(model_path)
            old, new = before.run_experiments(m, old_model), ME.run_experiments(m, model)
            save(f"{strategy}-{seed}-before.json", old); save(f"{strategy}-{seed}-after.json", new)
            assert old == new, (strategy, seed, "unexpected non-agent behavior change")
            summary["non_agent_pairs"].append({"strategy": strategy, "seed": seed, "exact_equal": True})
    if args.live_qwen:
        ctx = {"hypotheses": [{"id": "H_sec", "gene": "g6"}, {"id": "H_coup", "reactions": ["BIOMASS"]}],
               "evidence": [{"id": "E_sec", "type": "bgc", "hypothesis_id": "H_sec", "product_exchange": "EX_sec_e"},
                            {"id": "E_coup", "type": "bgc", "hypothesis_id": "H_coup", "product_exchange": "EX_coup_e"}]}
        m = manifest(ctx, strategy="qwen", seed=0, max_experiments=4, max_solver_calls=60)
        m["model"] = str(model_path)
        save("agent-manifest.json", m)
        for tag, module in (("before", before), ("after", ME)):
            chat = UsageChat(ME.resolve_config(m))
            result = module.run_experiments(m, old_model if tag == "before" else model, chat=chat)
            save(f"agent-{tag}.json", result)
            save(f"transcript-{tag}.json", chat.transcript)
            save(f"requests-{tag}.json", chat.log)
            summary["agent_runs"].append({"version": tag, "usage": chat.usage(),
                "selected": [x.get("experiment_id") for x in result["selection_log"] if x.get("experiment_id")],
                "fallback": any(x["source"] == "evidence_fallback" for x in result["selection_log"]),
                "valid_choices": all(x.get("valid") for x in result["selection_log"] if x["source"] == "qwen"),
                "solver_calls": result["solver"]["used"]})
            save("summary.json", summary)
            assert chat.usage()["ok"] > 0 and chat.usage()["failed"] == 0, "endpoint diagnostic failed"
            assert summary["agent_runs"][-1]["valid_choices"] and not summary["agent_runs"][-1]["fallback"]
            if tag == "after":
                for entry, transcript in zip(result["selection_log"], chat.transcript):
                    assert entry["prompt_sha256"] == ME.canonical_hash(transcript["messages"])
    save("summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
