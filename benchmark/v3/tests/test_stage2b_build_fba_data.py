import importlib.util
import json
import os
import pickle
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

V3_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = V3_ROOT.parents[1]
if str(V3_ROOT) not in sys.path:
    sys.path.insert(0, str(V3_ROOT))

from stage2b.build_fba_data import (
    build_fba_data_dir,
    parse_medium_csv_text,
    parse_seed_stoichiometry,
)


def _load_sidecar_module():
    sidecar_path = REPO_ROOT / "scripts" / "stoichiometric_locus_refiner.py"
    spec = importlib.util.spec_from_file_location("stoichiometric_locus_refiner", sidecar_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _write_synthetic_gapsmith_dir(root: Path) -> None:
    (root / "biomass").mkdir(parents=True, exist_ok=True)
    (root / "media").mkdir(parents=True, exist_ok=True)

    # Metabolites table
    (root / "seed_metabolites_edited.tsv").write_text(
        "id\tabbreviation\tname\tformula\tmass\tsource\tinchi\tcharge\tis_core\tis_obsolete\tlinked_compound\tis_cofactor\tdeltag\tdeltagerr\tpka\tpkb\tabstract_compound\tcomprised_of\taliases\tsmiles\tmnx_id\n"
        "cpd00001\th2o\tH2O\tH2O\t18\tPrimary\t\t0\t1\t0\t\t0\t0\t0\t\t\t\t\t\t\tMNXM2\n"
        "cpd00002\tatp\tATP\tC10H12N5O13P3\t503\tPrimary\t\t-4\t1\t0\t\t0\t0\t0\t\t\t\t\t\t\tMNXM3\n"
        "cpd00008\tadp\tADP\tC10H12N5O10P2\t424\tPrimary\t\t-3\t1\t0\t\t0\t0\t0\t\t\t\t\t\t\tMNXM7\n"
        "cpd00009\tpi\tPhosphate\tHO4P\t96\tPrimary\t\t-2\t1\t0\t\t0\t0\t0\t\t\t\t\t\t\tMNXM9\n"
        "cpd00020\tpyr\tPyruvate\tC3H3O3\t87\tPrimary\t\t-1\t1\t0\t\t0\t0\t0\t\t\t\t\t\t\tMNXM23\n"
        "cpd00027\tglc\tD-Glucose\tC6H12O6\t180\tPrimary\t\t0\t1\t0\t\t0\t0\t0\t\t\t\t\t\t\tMNXM41\n"
        "cpd00067\th\tH+\tH\t1\tPrimary\t\t1\t1\t0\t\t0\t0\t0\t\t\t\t\t\t\tMNXM1\n"
        "cpd99999\tdead\tDeadEndMet\tC3H3O3\t87\tPrimary\t\t-1\t1\t0\t\t0\t0\t0\t\t\t\t\t\t\tMNXM99999\n",
        encoding="utf-8",
    )

    # Reactions table:
    # rxn00001: D-Glucose transport (e -> c)
    # rxn00002: D-Glucose_c + ADP_c + Pi_c -> 2 Pyruvate_c + ATP_c + H2O_c
    # rxn00003: Dead-end branch Pyruvate_c -> DeadEndMet_c (no consumer of cpd99999_c)
    # rxn00004: Removed status reaction (should be filtered out)
    # rxn00005: Excluded bacterial reaction (should be filtered out)
    (root / "seed_reactions_corrected.tsv").write_text(
        "id\tabbreviation\tname\tcode\tstoichiometry\tis_transport\tequation\tdirection\treversibility\tstatus\tis_obsolete\tlinked_reaction\tnotes\tsource\tdeltag\tdeltagerr\tcompound_ids\tpathways\tec_numbers\tgapseq.status\n"
        'rxn00001\tglc_tr\tGlucose transport\t\t-1:cpd00027:1:0:"D-Glucose";1:cpd00027:0:0:"D-Glucose"\t1\t\t=\t=\tOK\t0\t\t\tPrimary\t0\t0\t\tTransport\t2.7.1.69\tapproved\n'
        'rxn00002\tglc_cat\tGlucose catabolism\t\t-1:cpd00027:0:0:"D-Glucose";2:cpd00020:0:0:"Pyruvate"\t0\t\t>\t>\tOK\t0\t\t\tPrimary\t0\t0\t\tGlycolysis\t2.7.1.1\tapproved\n'
        'rxn00003\tdead_rxn\tDead end branch\t\t-1:cpd00020:0:0:"Pyruvate";1:cpd99999:0:0:"DeadEndMet"\t0\t\t>\t>\tOK\t0\t\t\tPrimary\t0\t0\t\tOther\t1.1.1.999\tapproved\n'
        'rxn00004\trem_rxn\tRemoved reaction\t\t-1:cpd00020:0:0:"Pyruvate";1:cpd00027:0:0:"D-Glucose"\t0\t\t>\t>\tOK\t0\t\t\tPrimary\t0\t0\t\tOther\t1.1.1.998\tremoved\n'
        'rxn00005\texcl_rxn\tExcluded reaction\t\t-1:cpd00020:0:0:"Pyruvate";1:cpd00027:0:0:"D-Glucose"\t0\t\t>\t>\tOK\t0\t\t\tPrimary\t0\t0\t\tOther\t1.1.1.997\tapproved\n',
        encoding="utf-8",
    )

    (root / "biomass" / "excluded_reactions.tsv").write_text(
        "tax\trxn.id\trxn.name\treason\n"
        "Bacteria\trxn00005\tExcluded reaction\ttest exclusion\n",
        encoding="utf-8",
    )

    biomass_template = {
        "id": "Synthetic_Biomass",
        "name": "SyntheticBiomass",
        "domain": "Bacteria",
        "energy_GAM": 0.0,
        "met_groups": [
            {
                "group_name": "DNA",
                "mass": 0.1,
                "unit_components": "MOLFRACTION",
                "components": [{"id": "cpd00020", "comp": "c", "coef": 1.0, "link": ""}],
            },
            {
                "group_name": "RNA",
                "mass": 0.1,
                "unit_components": "MOLFRACTION",
                "components": [{"id": "cpd00020", "comp": "c", "coef": 1.0, "link": ""}],
            },
            {
                "group_name": "Protein",
                "mass": 0.3,
                "unit_components": "MOLFRACTION",
                "components": [{"id": "cpd00020", "comp": "c", "coef": 1.0, "link": ""}],
            },
        ],
    }
    (root / "biomass" / "biomass_Gram_neg.json").write_text(
        json.dumps(biomass_template), encoding="utf-8"
    )
    (root / "biomass" / "biomass_Gram_pos.json").write_text(
        json.dumps(biomass_template), encoding="utf-8"
    )

    (root / "seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv").write_text(
        "MS ID\tOld MS ID\tExternal ID\tSource\n"
        "rxn00001\t\t2.7.1.69\tEnzyme Class\n"
        "rxn00002\t\t2.7.1.1\tEnzyme Class\n"
        "rxn00003\t\t1.1.1.999\tEnzyme Class\n",
        encoding="utf-8",
    )
    (root / "altec.csv").write_text(
        "old_ec,new_ec\n"
        "2.7.1.1,2.7.1.2\n"
        "2.7.1.1,2.7.99.1\n",
        encoding="utf-8",
    )
    (root / "mnxref_seed.tsv").write_text(
        "rxn_mnx\tmets\tec\tseed\tbigg\tkegg\tmetacyc\n"
        "MNXR1\t\t2.7.1.69\trxn00001\t\t\t\n",
        encoding="utf-8",
    )
    (root / "media" / "MM_glu.csv").write_text(
        '"compounds","name","maxFlux"\n'
        '"cpd00027","D-Glucose, anhydrous",100\n'
        '"cpd00001","H2O",1000\n'
        '"cpd00009","Phosphate",100\n'
        '"cpd99999","ZeroFluxCompound",0\n',
        encoding="utf-8",
    )


class TestStage2BBuildFbaData(unittest.TestCase):
    def test_parse_seed_stoichiometry_and_medium_csv(self):
        raw = '-1:cpd00027:1:0:"D-Glucose, extracellular";1:cpd00027:0:0:"D-Glucose"'
        stoich, names = parse_seed_stoichiometry(raw)
        self.assertEqual(stoich, {"cpd00027_e": -1.0, "cpd00027_c": 1.0})
        self.assertEqual(names["cpd00027_e"], "D-Glucose, extracellular")

        csv_text = (
            '"compounds","name","maxFlux"\n'
            '"cpd00027","D-Glucose, alpha",100\n'
            '"cpd00001","H2O",0\n'
        )
        self.assertEqual(parse_medium_csv_text(csv_text), ["cpd00027"])

    def test_build_from_synthetic_dir_and_tar_and_run_sidecar(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref_dir = tmp_path / "gapsmith_ref"
            _write_synthetic_gapsmith_dir(ref_dir)

            out_dir = tmp_path / "fba_data_dir"
            receipt = build_fba_data_dir(
                gapsmith_source=ref_dir,
                output_dir=out_dir,
                reaction_filter="usable",
                tight_bounds_mode="fva",
                verify_pinned_sha256=False,
            )

            expected_files = (
                "universal_metadata.json",
                "fba_matrices_v6RE.pkl",
                "tight_bounds_v6_neg.pkl",
                "tight_bounds_v6_pos.pkl",
                "seedr2ec.pkl",
                "seedec2r.pkl",
                "medium.pkl",
                "receipt.json",
            )
            for fn in expected_files:
                self.assertTrue((out_dir / fn).is_file(), f"Missing {fn}")

            meta = json.loads((out_dir / "universal_metadata.json").read_text(encoding="utf-8"))
            allrxns = meta["allrxns"]
            # Internal usable non-excluded: rxn00001_c, rxn00002_c, rxn00003_c
            # Exchange: EX_cpd00027_e
            # Biomass: biomass_GmNeg, biomass_GmPos, EX_biomass
            self.assertIn("rxn00001_c", allrxns)
            self.assertIn("rxn00002_c", allrxns)
            self.assertIn("rxn00003_c", allrxns)
            self.assertNotIn("rxn00004_c", allrxns)
            self.assertNotIn("rxn00005_c", allrxns)
            self.assertIn("EX_cpd00027_e", allrxns)
            self.assertEqual(allrxns[-3:], ["biomass_GmNeg", "biomass_GmPos", "EX_biomass"])

            with open(out_dir / "tight_bounds_v6_neg.pkl", "rb") as f:
                tb_neg = pickle.load(f)
            idx_dead = allrxns.index("rxn00003_c")
            idx_active = allrxns.index("rxn00002_c")
            self.assertAlmostEqual(float(tb_neg["lb_tight"][idx_dead]), 0.0, places=7)
            self.assertAlmostEqual(float(tb_neg["ub_tight"][idx_dead]), 0.0, places=7)
            self.assertGreater(float(tb_neg["ub_tight"][idx_active]), 0.0)
            self.assertGreater(receipt["tight_bounds"]["neg"]["max_univ_bio_minimal"], 0.0)

            with open(out_dir / "seedec2r.pkl", "rb") as f:
                seedec2r = pickle.load(f)
            # altec expansion: 2.7.1.2 should inherit rxn00002 from 2.7.1.1,
            # while unspecific 2.7.99.1 should not pollute specific ECs.
            self.assertIn("rxn00002", seedec2r["2.7.1.1"])
            self.assertIn("rxn00002", seedec2r["2.7.1.2"])

            # Also test building from a tar archive with data_merged/ prefix
            tar_path = tmp_path / "data_merged.tar"
            with tarfile.open(tar_path, "w") as tf:
                for p in sorted(ref_dir.rglob("*")):
                    if p.is_file():
                        rel = p.relative_to(ref_dir)
                        tf.add(p, arcname=f"data_merged/{rel}")

            out_tar_dir = tmp_path / "fba_data_from_tar"
            receipt_tar = build_fba_data_dir(
                gapsmith_source=tar_path,
                output_dir=out_tar_dir,
                reaction_filter="usable",
                tight_bounds_mode="propagation",
                verify_pinned_sha256=False,
            )
            self.assertTrue(receipt_tar["gapsmith_source_is_tar"])
            self.assertEqual(
                receipt["output_sha256"]["universal_metadata.json"],
                receipt_tar["output_sha256"]["universal_metadata.json"],
            )

            # Prove the actual sidecar (scripts/stoichiometric_locus_refiner.py) loads and runs
            sidecar = _load_sidecar_module()
            refiner = sidecar.StoichiometricLocusRefiner(out_dir)
            agg = sidecar.SparseRoleAggregator()
            gene_ec_scores = {
                "P_GLC_TR": {"2.7.1.69": 0.90},
                "P_GLC_CAT": {"2.7.1.1": 0.25, "1.1.1.999": 0.40},
            }
            ec_probs, _ = agg.aggregate_ec_scores(gene_ec_scores)
            rxn_probs = agg.compute_reaction_probabilities(ec_probs, refiner.r2ec, refiner.rxn_meta)
            repair_res = refiner.solve_flux_certified_repair(
                rxn_probs, gram="neg", medium_name="minimal"
            )
            self.assertGreater(repair_res["max_univ_bio"], 0.0)
            self.assertGreater(repair_res["target_bio"], 0.0)
            self.assertIn("rxn00002_c", repair_res["repaired_rxns"])
            self.assertNotIn("rxn00003_c", repair_res["repaired_rxns"])
            updated = refiner.attribute_locus_responsibility(repair_res, gene_ec_scores, rxn_probs)
            # Active flux-carrying reaction EC 2.7.1.1 should be boosted above 0.25
            self.assertGreater(updated["P_GLC_CAT"]["2.7.1.1"], 0.25)

            sbml_path = tmp_path / "metabolic_model.xml"
            tsv_path = tmp_path / "reactions.tsv"
            refiner.export_model(repair_res, sbml_path, gene_ec_scores, "neg", "minimal")
            self.assertTrue(sbml_path.is_file() and sbml_path.stat().st_size > 100)
            self.assertTrue(tsv_path.is_file() and tsv_path.stat().st_size > 50)
            self.assertLessEqual(repair_res["model"]["max_mass_balance_residual"], 1e-6)


if __name__ == "__main__":
    unittest.main()
