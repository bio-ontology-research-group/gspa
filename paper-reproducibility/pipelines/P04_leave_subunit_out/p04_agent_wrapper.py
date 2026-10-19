#!/usr/bin/env python3
"""Run the Stage 3 agent sidecar with an explicit GO OBO for the symbolic guardrail.

Why: run_gspa_agent.py loads the GO DAG from REPO_ROOT/benchmark/neural/cafa6_recon/go-basic.obo
(REPO_ROOT = two directories above the sidecar).  On unimatrix01 the sidecar lives at
/data/hohndor/gspa/benchmark-unified/agent/, so it looked for
/data/hohndor/gspa/benchmark/neural/cafa6_recon/go-basic.obo, which does not exist.
SymbolicGuardrail then sets resources_loaded=False and every PROMOTE is rejected with
"Symbolic guardrail ontology/taxon resources are missing or incomplete; failing closed."
(the Oct-17 LSO run).  The taxon/anatomy tables themselves were present.

This wrapper imports the unmodified sidecar, points DEFAULT_GO_OBO at a real go.obo,
verifies that the guardrail loads (exits 3 otherwise, so a broken Stage 3 can no longer
pass silently), and then calls the sidecar's main() with the remaining arguments.

Usage: p04_agent_wrapper.py --go-obo /path/go.obo [--sidecar /path/run_gspa_agent.py] -- <sidecar args>
"""
import argparse
import importlib.util
import sys
from pathlib import Path

DEFAULT_SIDECAR = "/data/hohndor/gspa/benchmark-unified/agent/run_gspa_agent.py"


def main() -> int:
    argv = sys.argv[1:]
    if "--" in argv:
        i = argv.index("--")
        own, rest = argv[:i], argv[i + 1:]
    else:
        own, rest = argv, []
    ap = argparse.ArgumentParser()
    ap.add_argument("--go-obo", required=True, type=Path)
    ap.add_argument("--sidecar", default=DEFAULT_SIDECAR, type=Path)
    a = ap.parse_args(own)

    spec = importlib.util.spec_from_file_location("run_gspa_agent", str(a.sidecar))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.DEFAULT_GO_OBO = a.go_obo

    res_dir = None
    if "--resources-dir" in rest:
        res_dir = Path(rest[rest.index("--resources-dir") + 1])
    g = mod.SymbolicGuardrail(res_dir) if res_dir else mod.SymbolicGuardrail()
    if not g.resources_loaded:
        print(f"[p04-agent] guardrail resources failed to load (resources={res_dir}, obo={a.go_obo})",
              file=sys.stderr)
        return 3
    chk = g.check_go_term("GO:0000345", "NCBITaxon_2")
    print(f"[p04-agent] guardrail loaded: {len(g.known_taxa)} taxa, {len(g.known_go_terms)} GO terms; "
          f"self-test GO:0000345@NCBITaxon_2 valid={chk.get('valid')}", file=sys.stderr)

    sys.argv = [str(a.sidecar)] + rest
    return mod.main()


if __name__ == "__main__":
    sys.exit(main())
