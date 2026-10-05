# Flux constants of the cross-section bins: fills the 4th column of xs_real_bin.txt files (replaces og_calculate_num.C,
# get_flux_factor.C and get_flux_factor_3d.C).
#
#   python3 flux_constants.py <xs_real_bin.txt or directory> [...] [--pot 5e19] [--emax 8] [-v]
#
# A directory is searched recursively for xs_real_bin.txt files; each file is updated in place (-v prints every bin with
# the ranges used).
#
# constant = (integrated numu flux [numu/POT/cm2]) x POT x (number of Ar targets) x 1e-36 (pb) x (bin widths)
#   - the flux is integrated from 0 to --emax GeV (midpoint sum of TGraph::Eval of gh_averaged_numu_flux, as the old
#     macros);
#   - the bin widths are those of every differential variable of the bin (outer slice, inner slice, bin), energies and
#     momenta in GeV (cross sections in pb/GeV, pb/GeV^2, ...), cosines without units;
#   - neutrino energy bins: no width, the flux is integrated over the Enu range of the bin instead (the first / last bin
#     from 0 / up to --emax, where their events come from);
#   - bins without a variable (0p / Np, proton multiplicity, the 0p bin of the proton measurements): total flux only.
# The ranges are read from the bin names (get_xs_signal_no in cuts.h), e.g.
#   numuCC.Np.inside.Eavail.le.250.gt.150.costheta.le.0.5.gt.0.3.Emu.le.450.gt.350
# Open first / last bins (".le.e1" / ".gt.en"): extended by 2x the width of the neighbouring bin of the same slice (covers
# ~100% of the events of a first bin, 60-97% of a last bin), but not beyond the physical limits (Emu >= m_mu, Kp >= 45 MeV,
# energies and pt >= 0, cosines in [-1, 1]); without a neighbouring bin (a slice with one edge or none) the defaults
# below, from the 1D binnings.
import os, re, sys
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
FLUX_FILE, FLUX_GRAPH = os.path.join(HERE, 'gh_averaged_numu_flux.root'), 'gh_averaged_numu_flux'   # numu/POT/GeV/cm2
NSTEP = 1000

# number of argon targets in the fiducial volume
DENSITY, VOLUME, M_MOL, N_A = 1.3836, 5.82515e7, 39.95, 6.022e23   # g/cm3, cm3, g/mol, 1/mol
N_TARGET = DENSITY * VOLUME * N_A / M_MOL
PB = 1e-36   # cm2

# variable: (scale to the units of the cross section, physical lower limit, physical upper limit, default lower / upper
# limit of an open bin without a neighbouring bin) - names are in MeV
VARS = {'Emu': (1e-3, 105.66, None, 105.66, 2500), 'nu': (1e-3, 0, None, 0, 2500), 'Eavail': (1e-3, 0, None, 0, 1700),
        'Kp': (1e-3, 45, None, 45, 720), 'pt': (1e-3, 0, None, 0, 1050), 'pl': (1e-3, None, None, -450, 2700),
        'costheta': (1, -1, 1, -1, 1), 'costhetap': (1, -1, 1, -1, 1), 'costhetamup': (1, -1, 1, -1, 1)}
K_OPEN = 2   # open bins: at most this many widths of the neighbouring bin

NUM = r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?'
SEG = re.compile(r'([A-Za-z]+)((?:\.(?:le|gt)\.%s)*)(?:\.|$)' % NUM)

def segments(name):
    '''[(context, variable, le, gt)] of a bin name; le / gt None for an open side'''
    rest = name.split('.inside', 1)[1].lstrip('.')
    head = name[:len(name) - len(rest)]
    out, pos = [], 0
    while pos < len(rest):
        m = SEG.match(rest, pos)
        if not m: raise ValueError('cannot parse bin name ' + name)
        le = re.search(r'\.le\.(%s)' % NUM, m.group(2)); gt = re.search(r'\.gt\.(%s)' % NUM, m.group(2))
        out.append((head + rest[:pos], m.group(1), float(le.group(1)) if le else None, float(gt.group(1)) if gt else None))
        pos = m.end()
    return out

class Flux:
    def __init__(self):
        x, y = uproot.open(FLUX_FILE)[FLUX_GRAPH].values()
        self.x, self.y = np.asarray(x, float), np.asarray(y, float)
    def eval(self, e):
        '''TGraph::Eval: linear interpolation, linear extrapolation from the two closest points outside the graph'''
        i = np.clip(np.searchsorted(self.x, e) - 1, 0, len(self.x) - 2)
        x0, x1, y0, y1 = self.x[i], self.x[i + 1], self.y[i], self.y[i + 1]
        return y0 + (e - x0) * (y1 - y0) / (x1 - x0)
    def integral(self, lo, hi):
        e = lo + (hi - lo) * (np.arange(NSTEP) + 0.5) / NSTEP
        return self.eval(e).sum() * (hi - lo) / NSTEP

def open_range(var, le, gt, edges):
    '''[lo, hi] in MeV of a bin of var; edges: the sorted edges of its slice'''
    scale, pmin, pmax, dlo, dhi = VARS[var]
    if le is None and gt is None: return dlo, dhi   # a slice without bins: the whole range
    lo, hi = gt, le
    if lo is None:   # first bin
        k = edges.index(le)
        lo = le - K_OPEN * (edges[k + 1] - le) if k + 1 < len(edges) else dlo
        if pmin is not None: lo = max(lo, pmin)
        if pmax is not None: lo = pmin                # cosines: to -1
    if hi is None:   # last bin
        k = edges.index(gt)
        hi = gt + K_OPEN * (gt - edges[k - 1]) if k > 0 else dhi
        if pmax is not None: hi = pmax                # cosines: to 1
    if not hi > lo: raise ValueError('empty range for %s: %s %s' % (var, lo, hi))
    return lo, hi

def constants(names, flux, pot, emax):
    '''flux constant and the ranges used of every bin'''
    # edges of every slice: (context, variable) -> sorted edges
    slices = {}
    for nm in names:
        for ctx, var, le, gt in segments(nm):
            slices.setdefault((ctx, var), set()).update(v for v in (le, gt) if v is not None)
    slices = {k: sorted(v) for k, v in slices.items()}
    total = flux.integral(0, emax) * pot * N_TARGET * PB
    out = []
    for nm in names:
        c, ranges = total, []
        for ctx, var, le, gt in segments(nm):
            if var == 'Enu':   # flux over the Enu range of the bin instead of the width
                lo = gt if gt is not None else 0; hi = le if le is not None else emax * 1000
                c = c / flux.integral(0, emax) * flux.integral(lo / 1000, hi / 1000)
            else:
                if var not in VARS: raise ValueError('unknown variable %s in %s' % (var, nm))
                lo, hi = open_range(var, le, gt, slices[(ctx, var)])
                c *= (hi - lo) * VARS[var][0]
            ranges.append('%s [%g, %g]' % (var, lo, hi))
        out.append((c, ranges))
    return out

def update(fn, flux, pot, emax, verbose):
    lines = open(fn).read().split('\n')
    rows = [(i, l.split()) for i, l in enumerate(lines) if l.split() and not l.lstrip().startswith('#')]
    end = next((k for k, (i, t) in enumerate(rows) if t[0] == '-1'), len(rows))
    rows = rows[:end]
    if not rows or not all(re.match(r'numuCC\.(0p|Np|1p|2p|gt2p)\.inside(\.|$)', t[2]) for i, t in rows):
        print('%s: skipped, not get_xs_signal_no bin names' % fn); return
    res = constants([t[2] for i, t in rows], flux, pot, emax)
    for (i, t), (c, ranges) in zip(rows, res):
        t[3] = '%.6g' % c
        lines[i] = '%-6s  %-4s  %-90s  %-14s  %s  %s' % tuple(t[:6])
        if verbose: print('  %4s %-90s %12s  %s' % (t[0], t[2], t[3], ', '.join(ranges)))
    open(fn, 'w').write('\n'.join(lines))
    print('%s: %d bins' % (fn, len(rows)))

if __name__ == '__main__':
    args = sys.argv[1:]
    def opt(name, default):
        if name in args:
            k = args.index(name); v = args[k + 1]; del args[k:k + 2]; return float(v)
        return default
    pot, emax = opt('--pot', 5e19), opt('--emax', 8.)
    verbose = '-v' in args
    args = [a for a in args if a != '-v']
    files = []
    for a in args:
        if os.path.isdir(a):
            files += sorted(os.path.join(r, f) for r, _, fs in os.walk(a) for f in fs if f == 'xs_real_bin.txt')
        else: files.append(a)
    flux = Flux()
    print('POT %g, flux 0-%g GeV: %.6g numu/POT/cm2, %.6g Ar targets' % (pot, emax, flux.integral(0, emax), N_TARGET))
    for fn in files: update(fn, flux, pot, emax, verbose)
