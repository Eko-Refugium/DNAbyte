import sys
import os

# Ensure the repo root and tests directory are on sys.path so that
# imports like `from dnabyte...` and `from tests...` resolve correctly
# regardless of where pytest is invoked from.
repo_root = os.path.dirname(os.path.abspath(__file__))
tests_dir = os.path.join(repo_root, 'tests')
norec4dna_dir = os.path.join(repo_root, 'dnabyte', 'encoding', 'dna_aeon', 'DNA_Aeon', 'NOREC4DNA')

for path in (repo_root, tests_dir, norec4dna_dir):
    if path not in sys.path:
        sys.path.insert(0, path)
