"""Run offgrowth.runner while dumping all thread stacks to stderr every STACK_DUMP_SEC seconds (diagnostics only)."""
import faulthandler
import os
import sys

from offgrowth import runner

faulthandler.dump_traceback_later(int(os.environ.get("STACK_DUMP_SEC", "120")), repeat=True, file=sys.stderr)
sys.exit(runner.main())
