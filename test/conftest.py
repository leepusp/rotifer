"""
Shared pytest setup: make test/golden.py importable from every test
directory (``import golden``).
"""

import os
import sys

TEST = os.path.dirname(os.path.realpath(__file__))
if TEST not in sys.path:
    sys.path.insert(0, TEST)
