#!/usr/bin/env python3
# Entries of the checkout trees and the summed pot_tor875 of T_pot of ROOT files (local paths or xrootd URLs), as json
# on stdout: [[file, {tree: entries, ..., 'pot': ...}], ...] in the order given. The trees (and branches that must exist) are those every
# production checkout has (run_production.EVENT_TREES, T_pot, nuselection/SubRun); a missing tree or branch is an error.
# Needs PyROOT: run in the SL7 container with uboonecode set up, or on a grid worker (grid jobs write it next to their
# outputs, so the scripts that follow the jobs need no ROOT).
#   python3 file_info.py <file> [<file> ...] > info.json
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_production as rp

TREES = dict(rp.EVENT_TREES)
TREES[rp.POT_TREE[0]] = rp.POT_TREE[1]
TREES['nuselection/SubRun'] = ['run', 'subRun']


def main():
    out = [[p, rp.inspect_file(p, TREES, pot=True)] for p in sys.argv[1:]]
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
