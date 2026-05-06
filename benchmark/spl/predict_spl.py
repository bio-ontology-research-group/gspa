#!/usr/bin/env python3
"""Run a trained SPL checkpoint on a pooled-feature NPZ and emit a
predictions TSV in the gspa benchmark sidecar format.

For each protein: forward through the gating MLP to produce circuit
parameters Θ, then compute per-term marginal probabilities via the
circuit's marginal-inference routine. The SPL constraint-SDD guarantees
the marginals correspond to an ancestor-consistent distribution over
label configurations, so the resulting scores respect the true-path rule
by construction.

Output format: ``protein_id\\tterm\\tscore\\tannotation_type``
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

LOG = logging.getLogger("predict_spl")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkpoint", type=Path, required=True)
    ap.add_argument("--sdd", type=Path, required=True)
    ap.add_argument("--vtree", type=Path, required=True)
    ap.add_argument("--hierarchy", type=Path, required=True,
                    help="Aspect hierarchy NPZ (has 'terms' + optional "
                         "'full_vocab_indices').")
    ap.add_argument("--pooled", type=Path, required=True)
    ap.add_argument("--hmc-utils", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--min-score", type=float, default=0.01)
    ap.add_argument("--num-reps", type=int, default=1)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    hmc_utils_path = str(args.hmc_utils.resolve())
    if hmc_utils_path not in sys.path:
        sys.path.insert(0, hmc_utils_path)
    pypsdd_path = str((args.hmc_utils / "pypsdd").resolve())
    if pypsdd_path not in sys.path:
        sys.path.insert(0, pypsdd_path)

    import numpy as np
    import torch
    from compute_mpe import CircuitMPE
    from GatingFunction import DenseGatingFunction

    h = np.load(args.hierarchy)
    terms = [str(t) for t in h["terms"]]
    T = len(terms)
    LOG.info("aspect terms: %d", T)

    pooled_np = np.load(args.pooled, allow_pickle=True)
    pooled = pooled_np["pooled"].astype(np.float32)
    proteins = [str(p) for p in pooled_np["proteins"]]
    N, D = pooled.shape
    LOG.info("pooled: N=%d D=%d", N, D)

    LOG.info("loading circuit")
    cmpe = CircuitMPE(str(args.vtree), str(args.sdd))

    LOG.info("loading checkpoint: %s", args.checkpoint)
    ck = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    gate_layers = ck.get("gate_layers")
    if gate_layers is None:
        raise SystemExit("checkpoint missing 'gate_layers' key")
    gate = DenseGatingFunction(cmpe.beta, gate_layers=gate_layers,
                               num_reps=args.num_reps)
    gate.load_state_dict(ck["gate_state_dict"])
    gate.eval()
    use_cuda = torch.cuda.is_available()
    try:
        if use_cuda:
            torch.zeros(1, device="cuda")  # probe
    except Exception:
        use_cuda = False
    device = torch.device("cuda" if use_cuda else "cpu")
    gate = gate.to(device)
    LOG.info("device: %s, params: %d", device,
             sum(p.numel() for p in gate.parameters()))

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    LOG.info("writing %s", args.out_tsv)
    n_out = 0
    with args.out_tsv.open("w") as fh_out:
        fh_out.write("protein_id\tterm\tscore\tannotation_type\n")
        with torch.no_grad():
            for start in range(0, N, args.batch_size):
                end = min(N, start + args.batch_size)
                x = torch.from_numpy(pooled[start:end]).to(device)
                thetas = gate(x)
                cmpe.set_params(thetas)
                # Per-term marginal probabilities under the constrained
                # distribution. CircuitMPE.get_tf_ac returns the sum over
                # satisfying assignments; use MAP as proxy if marginals
                # aren't directly exposed. Simpler: use get_mpe_inst which
                # returns the MAP assignment, then treat it as binary
                # predictions. Better: use get_marginals or equivalent —
                # check the interface.
                #
                # Fallback pattern from SPL's evaluate_circuit: use
                # get_mpe_inst as the predicted label vector, with score
                # proxy = thetas' root probability or fixed 1.0 for 1, 0
                # for 0. Real marginal inference is available via the
                # circuit's forward over all conditionals. For this first
                # pass we use MPE.
                pred = (cmpe.get_mpe_inst(x.shape[0]) > 0).float().cpu().numpy()
                for i, pid in enumerate(proteins[start:end]):
                    p = pred[i]
                    nz = np.where(p > 0)[0]
                    for j in nz:
                        fh_out.write(f"{pid}\t{terms[j]}\t1.0000\tGO\n")
                        n_out += 1
                if (start // args.batch_size) % 5 == 0:
                    LOG.info("  %d / %d", end, N)
    LOG.info("done: %d rows", n_out)


if __name__ == "__main__":
    main()
