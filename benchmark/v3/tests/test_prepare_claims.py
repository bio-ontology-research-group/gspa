import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare_claims import deduplicate


class DeduplicationTests(unittest.TestCase):
    def test_only_identical_records_removed_preserving_distinct_homologs_and_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, out = Path(tmp) / "in.jsonl", Path(tmp) / "out.jsonl"
            claim = {"protein_id": "p", "function_type": "EC", "function_id": "1.1.1.1",
                     "source": "diamond", "raw_score": .8, "metadata": {"target": "homolog1"}}
            other_hit = {**claim, "metadata": {"target": "homolog2"}}
            other_score = {**claim, "raw_score": .9}
            rows = [claim, claim, other_hit, other_score, dict(reversed(list(claim.items())))]
            source.write_text("\n".join(json.dumps(row) for row in rows))
            before = source.read_bytes()
            audit = deduplicate(source, out)
            self.assertEqual(audit["identical_duplicates_removed"], 2)
            self.assertEqual([json.loads(l) for l in out.read_text().splitlines()], [claim, other_hit, other_score])
            self.assertEqual(source.read_bytes(), before)
            again = Path(tmp) / "again.jsonl"
            self.assertEqual(deduplicate(out, again)["identical_duplicates_removed"], 0)
            self.assertEqual(again.read_bytes(), out.read_bytes())

    def test_invalid_numeric_value_is_rejected_not_repaired(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, out = Path(tmp) / "in", Path(tmp) / "out"
            source.write_text('{"raw_score": NaN}\n')
            with self.assertRaises(ValueError):
                deduplicate(source, out)

    def test_existing_file_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, out = Path(tmp) / "in", Path(tmp) / "out"
            source.write_text('{}\n'); out.write_text('existing')
            with self.assertRaises(FileExistsError):
                deduplicate(source, out)
            self.assertEqual(out.read_text(), 'existing')
