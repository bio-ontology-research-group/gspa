import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_verified_arm as W  # noqa: E402


class RunVerifiedArm(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(dir=os.environ.get("TMPDIR")))

    def write(self, command):
        path = self.root / "command.json"
        path.write_text(json.dumps(command))
        return path

    def test_learned_flag_only_when_the_command_names_a_learned_model(self):
        for command, expected in ((["java", "-jar", "g.jar", "integrate"], False),
                                  (["java", "-jar", "g.jar", "integrate", "--learned-model", "m.json"], True)):
            path = self.write(command)
            with mock.patch.object(W.V, "run_batch", return_value={"returncode": 0, "artifacts": {},
                                                                   "spec": {"required_stages": []}}) as run:
                W.run(path, self.root / "manifest.tsv", self.root / "inference")
            args, kwargs = run.call_args
            self.assertEqual(args[0], command)
            self.assertEqual(kwargs["allow_unsupported_learned"], expected)
            self.assertIn(Path(W.__file__), kwargs["sources"])
            self.assertIn(path, kwargs["sources"])

    def test_non_argv_command_is_rejected(self):
        for bad in ({"cmd": "java"}, [], ["java", 1]):
            with self.assertRaises(ValueError):
                W.load_command(self.write(bad))


if __name__ == "__main__":
    unittest.main()
