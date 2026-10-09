# Compare two covariance output files (xf_cov_matrix XsFlux/cov_<n>.root, xs_cov_matrix XsFlux/cov_xs.root, ...)
# object by object: matrices (TMatrixD), vectors (TVectorD) and histograms (contents with under/overflow). Used to check
# that the files with the weights of all knobs (merge_xf ... all, #period 0 in xf_input.txt) give the same results as the
# one-file-per-knob setup.
#
#   python3 compare_cov_files.py <file A> <file B> [--tol 0]
#
# For each object: identical, or the largest absolute and relative difference; objects only in one file are listed.
# Exit code 0 if everything agrees within --tol (relative, default 0 = bit for bit), 1 otherwise.
import sys, argparse
import numpy as np
import uproot

parser = argparse.ArgumentParser()
parser.add_argument('file_a')
parser.add_argument('file_b')
parser.add_argument('--tol', type=float, default=0., help='allowed relative difference (default 0: identical)')
args = parser.parse_args()

def contents(f, key):
    '''values of a TMatrixD / TVectorD / TH1 / TH2 as a flat numpy array, None for other classes'''
    cls = f.classnames()[key]
    obj = f[key]
    if cls.startswith('TMatrixT') or cls.startswith('TVectorT'):
        return np.asarray(obj.member('fElements'), dtype=np.float64)
    if cls.startswith('TH1') or cls.startswith('TH2'):
        return np.asarray(obj.values(flow=True), dtype=np.float64).ravel()
    return None

fa, fb = uproot.open(args.file_a), uproot.open(args.file_b)
ka = {k.split(';')[0]: k for k in fa.classnames()}
kb = {k.split(';')[0]: k for k in fb.classnames()}
n_bad = 0
for name in sorted(set(ka) | set(kb)):
    if name not in ka or name not in kb:
        print('%-40s only in %s' % (name, 'A' if name in ka else 'B')); n_bad += 1; continue
    a, b = contents(fa, ka[name]), contents(fb, kb[name])
    if a is None or b is None:
        print('%-40s %s (not compared)' % (name, fa.classnames()[ka[name]])); continue
    if a.shape != b.shape:
        print('%-40s different sizes %s %s' % (name, a.shape, b.shape)); n_bad += 1; continue
    if np.array_equal(a, b):
        print('%-40s identical (%d values)' % (name, a.size)); continue
    diff = np.abs(a - b)
    scale = np.maximum(np.abs(a), np.abs(b))
    rel = np.where(scale > 0, diff / np.where(scale > 0, scale, 1), 0)
    ok = rel.max() <= args.tol
    n_bad += not ok
    print('%-40s %s: %d of %d values differ, max abs diff %.3g, max rel diff %.3g'
          % (name, 'within tol' if ok else 'DIFFERENT', (diff > 0).sum(), a.size, diff.max(), rel.max()))

print('\n%s' % ('all objects agree' if n_bad == 0 else '%d objects differ or are missing' % n_bad))
sys.exit(0 if n_bad == 0 else 1)
