#!/usr/bin/env python3
# Production of the BNB detector variation files from the unprocessed SURPRISE DetVar ntuples listed in
# numuCCana/mcc910wiki.txt, with the merge_det format that stores the CV once:
#
#   tree_trimmer <trimmer_config_det_pre.txt>      every CV and DetVar file: only the trees and branches bdt_convert and
#                                                  the analysis read (the branches of trimmer_config_numucc.txt, with
#                                                  pot_tor875), about 10 of the 79 GB of a 1M event file
#   bdt_convert -p1 -r1 -a<samdef> -t<config>      every pre-trimmed file
#   tree_trimmer -y<lo> -z<hi> -k0 -r0             split each file in run ranges of about --chunk-events CV events
#   merge_det <cv> <detvar> <out> -m2 -n<XX>       per run range: DetVar events matched to the CV, plus the flags of
#                                                  the CV events used; then hadd of the run ranges
#   hadd of the CV run ranges, merge_det <cv> <list> <out> -m3
#                                                  one CV file per CV sample, with flag_match_rse_XX of every DetVar file
#
# merge_det reads the DetVar file in the order of the CV events, which is very slow on whole files, so as the production
# team did it runs on run ranges (CV and DetVar events are matched within a run, and the POT ratios are per subrun, so
# the result is the one of the whole file). The run ranges of a CV sample come from the runs of its unprocessed file
# (cached in <state dir>/chunks.json); the CV ranges cover all runs, the DetVar ones the runs of the CV. The -m2 files
# and the -m3 CV file are in the order of the run ranges, the same for both.
#
# XX is the DetVar number of det_input.txt (#period; flag_match_rse_XX is what det_cov_matrix -rXX reads). The run 3
# DetVar files are merged twice, with the run 3a CV (output ..._3a) and with the run 3b CV (..._3b): the 1mil 3b CV for
# the ly* and WM* variations, the 500k one for recomb2 and sce. The bdt_convert and run range files are deleted once
# nothing needs them (the pre-trimmed file once its bdt_convert is done). At the end <state dir>/det_cv_input.txt lists every -m2 file with its CV file
# (configurations/det_cv_input.txt format).
#
# Checks: every pre-trimmed file has the entries of every tree it keeps and the POT of the unprocessed file; every run range file has the events and POT subruns of its source in that range (none lost or duplicated);
# every -m2 run range as for the whole file; the hadd files add up to their run ranges and to the bdt_convert CV file.
#
# Run inside the SL7 container with uboonecode set up, with bin/merge_det built from the det_disk_saver branch and
# bin/tree_trimmer with the POT of run ranges (-y -z) taken from the run limits:
#   python3 numuCCana/processing/run_production_det.py --dry-run
#   python3 numuCCana/processing/run_production_det.py -j 4
#   python3 numuCCana/processing/run_production_det.py --status
#   python3 numuCCana/processing/run_production_det.py --vars lya --runs 4d -j 2   # a subset (the -m3 files then only
#                                                                                   # have these flags, and are redone
#                                                                                   # with all of them on a full run)
#   python3 numuCCana/processing/run_production_det.py --reset 'm2:lya'          # forget the tasks matching a regex
#   (tasks of a run range are named <task>.<NN>, --status and --dry-run show them per task)
#
# Grid (on a gpvm, natively: no ROOT needed there; the apps come from the code tarball of --app-dir, which must have the
# TFile::Open versions of bdt_convert, tree_trimmer and merge_det; the BDT weights are shipped too):
#   python3 numuCCana/processing/run_production_det.py --grid --app-dir <build> [--wait]
#   python3 numuCCana/processing/run_production_det.py --grid --app-dir <build> --dry-run      # the first submission
# Every task is a grid job (onsite, 30 GB disk; memory / lifetime per kind in GRID_RES) that runs the command and then the
# task's check on the worker (--worker-check); each --grid run collects the finished jobs and submits the tasks whose
# dependencies are done (--wait: until everything is done). The first run submits jobs counting the events per run of
# the unprocessed CV files, for the run ranges (--worker-runs), and continues once they are in. Intermediate files go to
# /pnfs/uboone/scratch/users/<user>/grid/det/work (deleted once used), the final -m2 / -m3 files are written to
# .../grid/det/final and copied to --det-dir; job logs and markers in .../grid/det/jobs. A failed or lost job is
# resubmitted on the next run, up to --max-tries submissions; at most --max-detvar DetVar files have pre-trim / bdt_convert files at once.
#
# Resuming works as in run_production.py (whose helpers are used here): every task writes to a temporary file renamed
# after its check, its progress is in <state dir>/state.json, finished tasks are skipped unless their command (or -m3
# list) changed, and a deleted bdt_convert file is made again if a task that still has to run needs it.
import argparse, fcntl, json, os, re, shlex, shutil, signal, socket, subprocess, sys, time, traceback
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_production as rp
from run_production import StepError, State, log, tmp_name, remove, pid_alive, ROOT_LOCK
import grid as G

TOP = rp.TOP
OUT_DIR = '/pnfs/uboone/persistent/users/bbogart/numucc'

# CV samples: key, run, samdef (= unprocessed file name)
CVS = [
    dict(key='1', run='1', samdef='DetVar_Run123_v10_04_07_23_BNB_nu_overlay_cv_13a_surprise_reco2_hist_1'),
    dict(key='3a', run='3', samdef='DetVar_Run123_v10_04_07_23_BNB_nu_overlay_cv_13a_surprise_reco2_hist_3'),
    dict(key='3b_1mil', run='3', samdef='DetVar_Run123_v10_04_07_23_BNB_nu_overlay_cv_3b_1mil_surprise_reco2_hist'),
    dict(key='3b_500k', run='3', samdef='DetVar_Run123_v10_04_07_23_BNB_nu_overlay_cv_3b_500k_surprise_reco2_hist'),
    dict(key='4d', run='4d', samdef='DetVar_Run45_v10_04_07_19_BNB_nu_overlay_cv_surprise_reco2_hist_4d'),
    dict(key='5', run='5', samdef='DetVar_Run45_v10_04_07_19_BNB_nu_overlay_cv_surprise_reco2_hist_5'),
]
# variation name in the file names, DetVar number of det_input.txt
VARS = [('lyd', 1), ('lyr', 2), ('recomb2', 3), ('sce', 4), ('WMthetaXZ', 6), ('WMthetaYZ', 7), ('WMX', 8), ('WMYZ', 9), ('lya', 10)]
RUNS = ['1', '3', '4d', '5']
CHUNK_EVENTS = 40000  # CV events per run range (the production team's 50 runs of 4d were about 40k events)
RUN_MAX = 2147483647


def cv_keys(var, run):
    '''CV samples a DetVar file of this run is merged with, and the suffix of the -m2 output'''
    if run == '3':
        return [('3a', 'a'), ('3b_500k' if var in ('recomb2', 'sce') else '3b_1mil', 'b')]
    return [(run, '')]


def detvar_files(var, wiki):
    '''{run: unprocessed file} of a variation, from the wiki (run 1 / 3 files are named ..._reco2_1.root / _3.root)'''
    pat = re.compile(r'^DetVar_Run(123|45)_v[0-9_]+_BNB_nu_overlay_%s_surprise_reco2_(hist_)?(1|3|4d|5)\.root$' % var)
    out = {}
    for files in wiki:
        if len(files) != 1:
            continue
        m = pat.match(os.path.basename(files[0]))
        if m:
            if m.group(3) in out and out[m.group(3)] != files[0]:
                raise StepError('two unprocessed %s run %s files in the wiki: %s %s' % (var, m.group(3), out[m.group(3)], files[0]))
            out[m.group(3)] = files[0]
    return out


def samdef_of(path):
    '''samdef of a DetVar file: the file name, with _reco2_<run> written _reco2_hist_<run> as in the production samdefs'''
    return re.sub(r'_reco2_(1|3)$', r'_reco2_hist_\1', os.path.basename(path)[:-len('.root')])


def readable(path):
    '''the name ROOT opens a file by: on a grid worker (--worker-*), /pnfs files are read over xrootd'''
    if os.environ.get('DET_GRID_WORKER') and path.startswith('/pnfs/'):
        import grid_task
        return grid_task.xrootd_url(path)
    return path


def run_counts(path):
    '''{run: events} of wcpselection/T_eval'''
    with ROOT_LOCK:
        R = rp.root()
        f = R.TFile.Open(readable(path))
        if not f or f.IsZombie():
            raise StepError('can not open %s' % path)
        try:
            t = f.Get('wcpselection/T_eval')
            if not t:
                raise StepError('%s: no tree wcpselection/T_eval' % path)
            n = int(t.GetEntries())
            t.SetBranchStatus('*', 0)
            t.SetBranchStatus('run', 1)
            t.SetEstimate(n + 1)
            m = int(t.Draw('run', '', 'goff'))
            v = t.GetV1()
            out = {}
            for i in range(m):
                r = int(v[i])
                out[r] = out.get(r, 0) + 1
        finally:
            f.Close()
    if not out:
        raise StepError('%s: no events' % path)
    return out


def chunk_grid(counts, n):
    '''run ranges [lo, hi] of about n events (whole runs; a last range under n/2 events goes to the one before), the
    first starting at 0 and the last ending at RUN_MAX'''
    runs = sorted(counts)
    cuts = []
    k = 0
    for r in runs:
        k += counts[r]
        if k >= n:
            cuts.append(r)
            k = 0
    if k:
        if cuts and k < n / 2.:
            cuts[-1] = runs[-1]
        else:
            cuts.append(runs[-1])
    grid = []
    lo = 0
    for c in cuts:
        grid.append([lo, c])
        lo = c + 1
    grid[-1][1] = RUN_MAX
    return grid


class Grids:
    '''run ranges of the CV samples, from the runs of their unprocessed file, cached in <state dir>/chunks.json
    (fixed: take the cached ranges as they are, without looking at the files: grid workers)'''
    def __init__(self, path, n, fixed=False):
        self.path, self.n, self.fixed = path, n, fixed
        self.data = {}
        if os.path.exists(path):
            with open(path) as f:
                self.data = json.load(f)

    def cached(self, key, path):
        st = os.stat(path)
        c = self.data.get(key)
        if not c or c['input'] != path or c['size'] != st.st_size or c['mtime'] != int(st.st_mtime) or c['events'] != self.n:
            return None
        return c

    def store(self, key, path, counts):
        '''the run ranges of {run: events} of the unprocessed file path'''
        st = os.stat(path)
        c = dict(input=path, size=st.st_size, mtime=int(st.st_mtime), events=self.n, grid=chunk_grid(counts, self.n),
                 runs=[min(counts), max(counts)], total=sum(counts.values()))
        self.data[key] = c
        tmp = self.path + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(self.data, f, indent=1, sort_keys=True)
        os.replace(tmp, self.path)
        log('run ranges of %s: %d events in %d ranges' % (key, c['total'], len(c['grid'])))
        return c

    def get(self, key, path):
        if self.fixed:
            return self.data[key]
        c = self.cached(key, path)
        if not c:
            log('run ranges of %s: reading the runs of %s' % (key, path))
            c = self.store(key, path, run_counts(path))
        return c


# ---------------------------------------------------------------------------------------------------------------------
# ROOT queries
def query(path, trees, counts=None, sums=None, absent=()):
    '''{tree: entries} (each tree must exist with the branches listed), counts {key: (tree, selection)} -> entries passing,
    sums {key: (tree, expression, selection)} -> sum; trees in absent must not exist'''
    with ROOT_LOCK:
        R = rp.root()
        f = R.TFile.Open(readable(path))
        if not f or f.IsZombie():
            raise StepError('can not open %s' % path)
        info = {}
        try:
            for name, branches in trees.items():
                t = f.Get(name)
                if not t or not t.InheritsFrom('TTree'):
                    raise StepError('%s: no tree %s' % (path, name))
                missing = [b for b in branches if not t.GetBranch(b)]
                if missing:
                    raise StepError('%s: %s has no branch %s' % (path, name, ', '.join(missing)))
                info[name] = int(t.GetEntries())
            for key, (name, sel) in (counts or {}).items():
                info[key] = int(f.Get(name).GetEntries(sel))
            for key, (name, expr, sel) in (sums or {}).items():
                t = f.Get(name)
                n = int(t.GetEntries())
                total = 0.
                if n > 0:
                    t.SetEstimate(n + 1)
                    m = int(t.Draw(expr, sel, 'goff'))
                    v = t.GetV1()
                    total = sum(v[i] for i in range(m))
                info[key] = total
            for name in absent:
                if f.Get(name):
                    raise StepError('%s: unexpected tree %s' % (path, name))
        finally:
            f.Close()
        return info


def event_trees(suffix):
    return {'wcpselection/T_eval' + suffix: ['run', 'subrun', 'event', 'samdef', 'weight_cv', 'weight_spline', 'match_found'],
            'wcpselection/T_BDTvars' + suffix: ['numu_cc_flag', 'numu_score', 'all_veto_score', 'VtxAct_bdt_score'],
            'wcpselection/T_PFeval' + suffix: ['run', 'reco_nuvtxX', 'mcs_emu_MCS'],
            'wcpselection/T_KINEvars' + suffix: ['kine_reco_Enu'],
            'wcpselection/T_spacepoints' + suffix: ['Trecchargeblob_spacepoints_q'],
            'nuselection/NeutrinoSelectionFilter' + suffix: ['run', 'sub', 'evt'],
            'lantern/EventTree' + suffix: ['run', 'subrun', 'event', 'haveReco']}


def same_entries(info, trees, path, empty=False):
    n = [info[t] for t in trees]
    if n[0] <= 0 and not empty:
        raise StepError('%s: no events' % path)
    if len(set(n)) != 1:
        raise StepError('%s: event trees with different entries %s' % (path, dict((t, info[t]) for t in trees)))


def check_bdt(path, done):
    trees = event_trees('')
    info = query(path, dict(trees, **{'wcpselection/T_pot': ['runNo', 'subRunNo', 'pot_tor875']}),
                 sums={'pot': ('wcpselection/T_pot', 'pot_tor875', '')})
    same_entries(info, trees, path)
    if info['wcpselection/T_pot'] <= 0 or info['pot'] <= 0:
        raise StepError('%s: POT %g in %d subruns' % (path, info['pot'], info['wcpselection/T_pot']))
    return info


PRE_TREES = {'wcpselection/T_eval': ['run', 'subrun', 'event', 'weight_cv', 'match_found'],
             'wcpselection/T_BDTvars': ['numu_cc_flag'],
             'wcpselection/T_PFeval': ['run', 'reco_nuvtxX'],
             'wcpselection/T_KINEvars': ['kine_reco_Enu'],
             'wcpselection/T_spacepoints': ['Trecchargeblob_spacepoints_q', 'Trecchargeblob_spacepoints_real_cluster_id'],
             'nuselection/NeutrinoSelectionFilter': ['run', 'sub', 'evt'],
             'lantern/EventTree': ['run', 'subrun', 'event', 'haveReco'],
             'singlephotonana/eventweight_tree': ['GTruth_ResNum'],
             'wcpselection/T_pot': ['runNo', 'subRunNo', 'pot_tor875'],
             'nuselection/SubRun': ['run', 'subRun'],
             'lantern/potTree': ['run', 'subrun'],
             'singlephotonana/run_subrun_tree': ['run', 'subrun']}


def check_pre(src):
    '''pre-trim of the unprocessed file src: every tree kept has all its entries, the same POT, the big trees are gone'''
    def check(path, done):
        sums = {'pot': ('wcpselection/T_pot', 'pot_tor875', '')}
        info = query(path, PRE_TREES, sums=sums, absent=['singlephotonana/vertex_tree', 'shrreco3d/_dedx_tree'])
        same_entries(info, [t for t in event_trees('') if t in PRE_TREES], path)
        ref = query(src, dict((t, []) for t in PRE_TREES), sums=sums)
        bad = ['%s %d (%d)' % (t, info[t], ref[t]) for t in PRE_TREES if info[t] != ref[t]]
        if bad or not close(info['pot'], ref['pot']):
            raise StepError('%s: entries (unprocessed file) %s, POT %g (%g)' % (path, ', '.join(bad), info['pot'], ref['pot']))
        return info
    return check


def close(x, y):
    return abs(x - y) <= 1e-9 * max(abs(x), abs(y))


def check_trim(src, lo, hi, cv):
    '''run range [lo, hi] of src: the same events and POT subruns as src in that range (a DetVar range can be empty)'''
    def check(path, done):
        trees = event_trees('')
        sel, psel = 'run>=%d && run<=%d' % (lo, hi), 'runNo>=%d && runNo<=%d' % (lo, hi)
        info = query(path, dict(trees, **{'wcpselection/T_pot': ['runNo', 'subRunNo', 'pot_tor875']}),
                     counts={'n_out': ('wcpselection/T_eval', '!(%s)' % sel), 'n_pot_out': ('wcpselection/T_pot', '!(%s)' % psel)},
                     sums={'pot': ('wcpselection/T_pot', 'pot_tor875', '')})
        same_entries(info, trees, path, empty=not cv)
        ref = query(src, {}, counts={'n': ('wcpselection/T_eval', sel), 'n_pot': ('wcpselection/T_pot', psel)},
                    sums={'pot': ('wcpselection/T_pot', 'pot_tor875', psel)})
        if info['n_out'] or info['n_pot_out']:
            raise StepError('%s: %d events / %d subruns outside the runs %d-%d' % (path, info['n_out'], info['n_pot_out'], lo, hi))
        if info['wcpselection/T_eval'] != ref['n']:
            raise StepError('%s: %d events, %d in the runs %d-%d of %s' % (path, info['wcpselection/T_eval'], ref['n'], lo, hi, src))
        if info['wcpselection/T_pot'] != ref['n_pot'] or not close(info['pot'], ref['pot']):
            raise StepError('%s: POT %g in %d subruns, %g in %d in the runs %d-%d of %s (is bin/tree_trimmer built with the '
                            'POT of -y -z from the run limits?)' % (path, info['pot'], info['wcpselection/T_pot'], ref['pot'],
                                                                     ref['n_pot'], lo, hi, src))
        return info
    return check


def check_sum(info, parts, keys, done, path):
    '''a hadd file adds up to its parts'''
    for k in keys:
        x = sum(done[t]['info'][k] for t in parts)
        if not (close(info[k], x) if isinstance(x, float) else info[k] == x):
            raise StepError('%s: %s %s, the run ranges add up to %s' % (path, k, info[k], x))


def check_cv(bdt_task, parts):
    '''hadd of the CV run ranges: all the events and POT of the bdt_convert CV file'''
    def check(path, done):
        trees = event_trees('')
        info = query(path, dict(trees, **{'wcpselection/T_pot': ['runNo', 'subRunNo', 'pot_tor875']}),
                     sums={'pot': ('wcpselection/T_pot', 'pot_tor875', '')})
        same_entries(info, trees, path)
        check_sum(info, parts, list(trees) + ['wcpselection/T_pot', 'pot'], done, path)
        bdt = done[bdt_task]['info']
        if info['wcpselection/T_eval'] != bdt['wcpselection/T_eval'] or info['wcpselection/T_pot'] != bdt['wcpselection/T_pot'] \
                or not close(info['pot'], bdt['pot']):
            raise StepError('%s: %d events, POT %g in %d subruns; the bdt_convert file has %d, %g in %d' % (
                path, info['wcpselection/T_eval'], info['pot'], info['wcpselection/T_pot'], bdt['wcpselection/T_eval'], bdt['pot'],
                bdt['wcpselection/T_pot']))
        return info
    return check


def check_m2(det_no, cv_task, parts=None, empty=False):
    '''cv_task: the CV file merged (for the hadd of the run ranges, the bdt_convert CV file), parts: the run ranges hadd'ed,
    empty: a run range without DetVar events is fine'''
    def check(path, done):
        cv = done[cv_task]['info']
        trees = event_trees('_det')
        info = query(path, dict(trees, **{'wcpselection/T_pot_det': ['runNo'],
                                          'wcpselection/T_cv_rse_flag': ['run', 'subrun', 'event', 'flag', 'det_no'],
                                          'wcpselection/T_cv_pot_flag': ['runNo', 'subRunNo', 'flag', 'pot_tor875']}),
                     counts={'n_flag': ('wcpselection/T_cv_rse_flag', 'flag==1'),
                             'n_bad_det_no': ('wcpselection/T_cv_rse_flag', 'det_no!=%d' % det_no),
                             'n_pot_flag': ('wcpselection/T_cv_pot_flag', 'flag==1')},
                     sums={'pot': ('wcpselection/T_cv_pot_flag', 'pot_tor875', 'flag==1')},
                     absent=['wcpselection/T_eval_cv'])
        same_entries(info, trees, path, empty)
        if parts:
            check_sum(info, parts, list(trees) + ['wcpselection/T_pot_det', 'wcpselection/T_cv_rse_flag', 'wcpselection/T_cv_pot_flag',
                                                  'n_flag', 'n_pot_flag', 'pot'], done, path)
        if info['wcpselection/T_cv_rse_flag'] != cv['wcpselection/T_eval'] or info['wcpselection/T_cv_pot_flag'] != cv['wcpselection/T_pot']:
            raise StepError('%s: flags of %d events / %d subruns, the CV file has %d / %d' % (
                path, info['wcpselection/T_cv_rse_flag'], info['wcpselection/T_cv_pot_flag'], cv['wcpselection/T_eval'], cv['wcpselection/T_pot']))
        if info['n_flag'] != info['wcpselection/T_eval_det']:
            raise StepError('%s: %d CV events flagged, %d DetVar events' % (path, info['n_flag'], info['wcpselection/T_eval_det']))
        if info['n_bad_det_no']:
            raise StepError('%s: det_no is not %d' % (path, det_no))
        if (info['n_pot_flag'] <= 0 or info['pot'] <= 0) and not (empty and info['n_flag'] == 0):
            raise StepError('%s: POT %g in %d subruns' % (path, info['pot'], info['n_pot_flag']))
        return info
    return check


def check_m3(cv_task, m2s):
    '''m2s: [(det_no, m2 task)]'''
    def check(path, done):
        cv = done[cv_task]['info']
        trees = event_trees('_cv')
        nos = [n for n, _ in m2s]
        trees['wcpselection/T_eval_cv'] = trees['wcpselection/T_eval_cv'] + ['flag_match_rse_%d' % n for n in nos]
        pot_tree = {'wcpselection/T_pot_cv': ['runNo', 'subRunNo', 'pot_tor875'] + ['flag_match_rse_%d' % n for n in nos] + ['pot_tor875_rse_%d' % n for n in nos]}
        counts = dict(('n_flag_%d' % n, ('wcpselection/T_eval_cv', 'flag_match_rse_%d==1' % n)) for n in nos)
        sums = dict(('pot_%d' % n, ('wcpselection/T_pot_cv', 'pot_tor875_rse_%d' % n, 'flag_match_rse_%d==1' % n)) for n in nos)
        info = query(path, dict(trees, **pot_tree), counts=counts, sums=sums)
        same_entries(info, trees, path)
        if info['wcpselection/T_eval_cv'] != cv['wcpselection/T_eval'] or info['wcpselection/T_pot_cv'] != cv['wcpselection/T_pot']:
            raise StepError('%s: %d events / %d subruns, the CV file has %d / %d' % (
                path, info['wcpselection/T_eval_cv'], info['wcpselection/T_pot_cv'], cv['wcpselection/T_eval'], cv['wcpselection/T_pot']))
        for n, t in m2s:
            m2 = done[t]['info']
            if info['n_flag_%d' % n] != m2['n_flag']:
                raise StepError('%s: %d events flagged for %d, %d in %s' % (path, info['n_flag_%d' % n], n, m2['n_flag'], t))
            if not close(info['pot_%d' % n], m2['pot']):
                raise StepError('%s: POT %g for %d, %g in %s' % (path, info['pot_%d' % n], n, m2['pot'], t))
        return info
    return check


# ---------------------------------------------------------------------------------------------------------------------
# tasks: name -> dict(cmd, inputs, outputs, deps, check, content (file written before the command), intermediate)
def build_tasks(wiki, a, get_grid):
    tasks = {}
    order = []

    def add(name, **kw):
        kw.setdefault('content', None)
        kw.setdefault('intermediate', False)
        kw['name'] = name
        tasks[name] = kw
        order.append(name)

    app = lambda x: os.path.join(a.bin_dir, x)
    flags = lambda samdef: ['-p1', '-r1', '-a' + samdef, '-t' + a.config]
    bdt_out = lambda samdef: os.path.join(a.work_dir, 'bdt_convert_%s.root' % samdef)
    pre_out = lambda samdef: os.path.join(a.work_dir, 'pre_%s.root' % samdef)
    pre = lambda src: [app('tree_trimmer'), src, '{tmp}', a.pre_config, '-o1', '-k0', '-r0', '-v100000']
    chunk_out = lambda name, lo, hi: os.path.join(a.work_dir, 'ranges', '%s_%d_%d.root' % (name, lo, hi))
    trim = lambda src, lo, hi: [app('tree_trimmer'), src, '{tmp}', a.trim_config, '-o1', '-k0', '-r0', '-y%d' % lo, '-z%d' % hi, '-v100000']

    vars_sel = [(v, n) for v, n in VARS if not a.vars or v in a.vars.split(',')]
    runs_sel = [r for r in RUNS if not a.runs or r in a.runs.split(',')]
    cv_by_key = dict((c['key'], c) for c in CVS)
    wiki_cv = dict((os.path.basename(f[0])[:-len('.root')], f[0]) for f in wiki if len(f) == 1)

    # DetVar files and the CV samples they need
    m2_of_cv = {}
    det_jobs = []
    for var, no in vars_sel:
        files = detvar_files(var, wiki)
        for run in runs_sel:
            if run not in files:
                continue
            for key, suffix in cv_keys(var, run):
                m2_of_cv.setdefault(key, [])
            det_jobs.append((var, no, run, files[run]))

    # CV: bdt_convert, run ranges
    grids = {}
    for key in [c['key'] for c in CVS if c['key'] in m2_of_cv]:
        c = cv_by_key[key]
        if c['samdef'] not in wiki_cv:
            raise StepError('no unprocessed %s in the wiki' % c['samdef'])
        out = bdt_out(c['samdef'])
        cvt = 'bdt:cv_' + key
        src, pout = wiki_cv[c['samdef']], pre_out(c['samdef'])
        add('pre:cv_' + key, cmd=pre(src), inputs=[src], outputs=[pout], deps=[], check=check_pre(src), intermediate=True)
        add(cvt, cmd=[app('bdt_convert'), pout, '{tmp}'] + flags(c['samdef']),
            inputs=[pout], outputs=[out], deps=['pre:cv_' + key], check=check_bdt, intermediate=True)
        grids[key] = get_grid(key, wiki_cv[c['samdef']])
        for i, (lo, hi) in enumerate(grids[key]['grid']):
            add('trim:cv_%s.%02d' % (key, i), cmd=trim(out, lo, hi), inputs=[out],
                outputs=[chunk_out('bdt_convert_' + c['samdef'], lo, hi)], deps=[cvt], check=check_trim(out, lo, hi, True),
                intermediate=True)

    # DetVar: bdt_convert, then per CV sample: run ranges (those of the CV, within its runs), -m2 per range, hadd
    for var, no, run, f in det_jobs:
        samdef = samdef_of(f)
        out = bdt_out(samdef)
        bt = 'bdt:%s_%s' % (var, run)
        pt = 'pre:%s_%s' % (var, run)
        add(pt, cmd=pre(f), inputs=[f], outputs=[pre_out(samdef)], deps=[], check=check_pre(f), intermediate=True)
        add(bt, cmd=[app('bdt_convert'), pre_out(samdef), '{tmp}'] + flags(samdef), inputs=[pre_out(samdef)], outputs=[out],
            deps=[pt], check=check_bdt, intermediate=True)
        for key, suffix in cv_keys(var, run):
            rmin, rmax = grids[key]['runs']
            name = 'm2:%s_%s%s' % (var, run, suffix)
            parts = []
            for i, (lo, hi) in enumerate(grids[key]['grid']):
                dlo, dhi = max(lo, rmin), min(hi, rmax)
                ct = 'trim:cv_%s.%02d' % (key, i)
                tt = 'trim:%s_%s%s.%02d' % (var, run, suffix, i)
                add(tt, cmd=trim(out, dlo, dhi), inputs=[out], outputs=[chunk_out('bdt_convert_%s_%s' % (samdef, key), dlo, dhi)],
                    deps=[bt], check=check_trim(out, dlo, dhi, False), intermediate=True)
                mt = '%s.%02d' % (name, i)
                add(mt, cmd=[app('merge_det'), tasks[ct]['outputs'][0], tasks[tt]['outputs'][0], '{tmp}', '-m2', '-n%d' % no, '-t' + a.config],
                    inputs=[tasks[ct]['outputs'][0], tasks[tt]['outputs'][0]],
                    outputs=[chunk_out('merge_det_%s%s' % (samdef, suffix), lo, hi)], deps=[ct, tt], check=check_m2(no, ct, empty=True),
                    intermediate=True)
                parts.append(mt)
            files = [tasks[t]['outputs'][0] for t in parts]
            add(name, cmd=['hadd', '-f', '{tmp}'] + files, inputs=files,
                outputs=[os.path.join(a.det_dir, 'merge_det_%s%s.root' % (samdef, suffix))], deps=parts,
                check=check_m2(no, 'bdt:cv_' + key, parts=parts))
            m2_of_cv[key].append((no, name))

    # CV: hadd of the run ranges (in the order of the -m2 files; made once these are done), -m3
    for key in [c['key'] for c in CVS if c['key'] in m2_of_cv]:
        c = cv_by_key[key]
        m2s = sorted(m2_of_cv[key])
        parts = ['trim:cv_%s.%02d' % (key, i) for i in range(len(grids[key]['grid']))]
        files = [tasks[t]['outputs'][0] for t in parts]
        cvt = 'cv:' + key
        add(cvt, cmd=['hadd', '-f', '{tmp}'] + files, inputs=files,
            outputs=[os.path.join(a.work_dir, 'bdt_convert_%s_ranges.root' % c['samdef'])], deps=parts + [t for _, t in m2s],
            check=check_cv('bdt:cv_' + key, parts), intermediate=True)
        lst = os.path.join(a.work_dir, 'det_list_%s.txt' % c['samdef'])
        content = ''.join('%d %s\n' % (n, tasks[t]['outputs'][0]) for n, t in m2s) + 'end end\n'
        add('m3:cv_' + key, cmd=[app('merge_det'), tasks[cvt]['outputs'][0], lst, '{tmp}', '-m3', '-t' + a.config],
            inputs=[tasks[cvt]['outputs'][0]] + [tasks[t]['outputs'][0] for _, t in m2s],
            outputs=[os.path.join(a.det_dir, 'merge_det_cv_%s.root' % c['samdef'])], deps=[cvt] + [t for _, t in m2s],
            check=check_m3(cvt, m2s), content=(lst, content))
    return [tasks[n] for n in order]


def signature(t):
    s = ' '.join(x.replace('{tmp}', tmp_name(t['outputs'][0])) for x in t['cmd'])
    if t['content']:
        s += ' | ' + t['content'][1]
    return s


# ---------------------------------------------------------------------------------------------------------------------
def run_task(t, a, state):
    name = t['name']
    st = state.data.setdefault(name, {})
    logfile = os.path.join(a.state_dir, 'logs', re.sub(r'[^A-Za-z0-9_.-]', '_', name) + '.log')
    t0 = time.time()
    with state.lock:
        st.update(status='running', cmd=signature(t), started=time.strftime('%Y-%m-%d %H:%M:%S'), host=socket.gethostname(),
                  pid=None, log=logfile, error=None, deleted=False)
        state.save()
    log('%s: start' % name)
    try:
        for o in t['outputs']:
            os.makedirs(os.path.dirname(o), exist_ok=True)
            remove(tmp_name(o))
            remove(o)
        if t['content']:
            remove(t['content'][0])
            with open(t['content'][0], 'w') as f:
                f.write(t['content'][1])
        tmp = tmp_name(t['outputs'][0])
        cmd = [c.replace('{tmp}', tmp) for c in t['cmd']]
        with open(logfile, 'w') as lf:
            lf.write(' '.join(cmd) + '\n\n')
            lf.flush()
            p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=a.run_dir, start_new_session=True)
            rp.RUNNING[p.pid] = p
            with state.lock:
                st['pid'] = p.pid
                state.save()
            rc = p.wait()
            rp.RUNNING.pop(p.pid, None)
        if rp.STOP.is_set():
            raise StepError('interrupted')
        if rc != 0:
            raise StepError('exit code %d, see %s' % (rc, logfile))
        if not os.path.exists(tmp) or os.path.getsize(tmp) == 0:
            raise StepError('no output %s, see %s' % (tmp, logfile))
        info = t['check'](tmp, state.data)
        os.replace(tmp, t['outputs'][0])
        if t['content']:
            remove(t['content'][0])
        with state.lock:
            st.update(status='done', pid=None, info=info, seconds=round(time.time() - t0), finished=time.strftime('%Y-%m-%d %H:%M:%S'))
            state.save()
        log('%s: done (%d s)' % (name, time.time() - t0))
        return True
    except Exception as e:
        msg = str(e) if isinstance(e, StepError) else traceback.format_exc()
        with state.lock:
            st.update(status='interrupted' if rp.STOP.is_set() else 'failed', pid=None, error=msg)
            state.save()
        log('%s: %s: %s' % (name, 'stopped' if rp.STOP.is_set() else 'FAILED', msg.strip()))
        return False


def detvar_on_disk(tasks, state, todo, consumers, busy):
    '''DetVar files (not the CV ones) with a pre-trim or bdt_convert file being made (busy(name)), or made and still
    needed; at most --max-detvar at once (a new one starts with its pre-trim)'''
    def holds(t):
        n = t['name']
        return busy(n) or (n not in todo and not state.data.get(n, {}).get('deleted') and
                           any(u['name'] in todo for u in consumers.get(t['outputs'][0], [])))
    by_file = {}
    for t in tasks:
        if t['name'].startswith(('pre:', 'bdt:')) and not t['name'].startswith(('pre:cv_', 'bdt:cv_')):
            by_file.setdefault(t['name'].split(':', 1)[1], []).append(t)
    return sum(1 for ts in by_file.values() if any(holds(t) for t in ts))


def starts_detvar(t):
    return t['name'].startswith('pre:') and not t['name'].startswith('pre:cv_')


def delete_unneeded(tasks, a, state, todo, consumers):
    '''delete the bdt_convert and run range files whose consumers are all done (and not to be run again)'''
    if a.keep_work:
        return
    for p in tasks:
        if not p['intermediate'] or state.data.get(p['name'], {}).get('deleted'):
            continue
        f = p['outputs'][0]
        users = consumers.get(f, [])
        if users and all(t['name'] not in todo and state.data.get(t['name'], {}).get('status') == 'done' for t in users) and os.path.exists(f):
            remove(f)
            with state.lock:
                state.data.setdefault(p['name'], {})['deleted'] = True
                state.save()
            log('deleted %s' % f)


# ---------------------------------------------------------------------------------------------------------------------
# grid (--grid): every task is a grid job; its check runs on the worker (which has ROOT) right after the command
# (--worker-check: the same task list, built from the shipped wiki and run ranges), so the machine following the jobs
# needs no ROOT. Intermediate files (bdt_convert, run ranges, CV hadd) live in --work-dir on /pnfs scratch, the final
# -m2 / -m3 files are written to scratch and copied to --det-dir. Run ranges: the run counts of the unprocessed CV
# files come from grid jobs too (--worker-runs) the first time.
GRID_RES = {'runs': ('2GB', '2h'), 'pre': ('8GB', '16h'), 'bdt': ('4GB', '24h'), 'trim': ('6GB', '4h'), 'm2r': ('2GB', '4h'), 'hadd': ('2GB', '8h'),
            'm3': ('4GB', '8h')}


def grid_kind(name):
    if name.startswith('pre:'):
        return 'pre'
    if name.startswith('bdt:'):
        return 'bdt'
    if name.startswith('trim:'):
        return 'trim'
    if name.startswith('m2:'):
        return 'm2r' if re.search(r'\.\d+$', name) else 'hadd'
    if name.startswith('cv:'):
        return 'hadd'
    return 'm3'


def safe_name(name):
    return re.sub(r'[^A-Za-z0-9_.-]', '_', name)


def needed_cvs(wiki, a):
    '''[(CV key, unprocessed CV file)] of the DetVar files selected'''
    wiki_cv = dict((os.path.basename(f[0])[:-len('.root')], f[0]) for f in wiki if len(f) == 1)
    keys = set()
    for var, no in VARS:
        if a.vars and var not in a.vars.split(','):
            continue
        files = detvar_files(var, wiki)
        for run in RUNS:
            if (not a.runs or run in a.runs.split(',')) and run in files:
                keys.update(k for k, _ in cv_keys(var, run))
    out = []
    for c in CVS:
        if c['key'] in keys:
            if c['samdef'] not in wiki_cv:
                raise StepError('no unprocessed %s in the wiki' % c['samdef'])
            out.append((c['key'], wiki_cv[c['samdef']]))
    return out


def worker_args(a):
    '''the options a worker needs to build the same task list (paths of the shipped wiki and run ranges)'''
    args = ['--work-dir', a.work_dir, '--det-dir', a.det_dir, '--wiki', '{code}/processing/mcc910wiki.txt',
            '--pre-config', '{code}/processing/' + os.path.basename(a.pre_config),
            '--chunks', '{code}/processing/det_chunks.json', '--chunk-events', str(a.chunk_events)]
    if a.vars:
        args += ['--vars', a.vars]
    if a.runs:
        args += ['--runs', a.runs]
    return args


def grid_task_of(t, a, by_name, state, stage_dir, info_dir):
    '''the grid.py task of a task: its command with the apps, configs and inputs of the worker, then --worker-check'''
    sub = {}
    for x in ('bdt_convert', 'tree_trimmer', 'merge_det'):
        sub[os.path.join(a.bin_dir, x)] = '{code}/bin/' + x
    sub[a.trim_config] = '{code}/processing/' + os.path.basename(a.trim_config)
    sub[a.pre_config] = '{code}/processing/' + os.path.basename(a.pre_config)
    sub['-t' + a.config] = '-t{code}/processing/' + os.path.basename(a.config)
    sub['{tmp}'] = '{out0}'
    files = {}
    if t['content']:   # the -m3 list: the -m2 files read over xrootd
        import grid_task
        sub[t['content'][0]] = 'det_list.txt'
        files['det_list.txt'] = ''.join('%s %s\n' % (l.split()[0], grid_task.xrootd_url(l.split()[1])) if l.split()[0] != 'end' else l + '\n'
                                        for l in t['content'][1].strip().split('\n'))
    inputs = []
    for x in t['cmd']:
        if x in t['inputs'] and x not in sub:
            sub[x] = '{in%d}' % len(inputs)
            inputs.append(dict(path=x, mode='copy' if grid_kind(t['name']) == 'm2r' else 'stream'))
    cmd = [sub.get(x, x) for x in t['cmd']]
    if grid_kind(t['name']) == 'pre':
        # the raw files streamed over xrootd: without the TTreeCache, whose vector reads of the large baskets of the run 1 / 3
        # files dCache refuses ("[3002] Single readv transfer is too large", then garbage reads and runaway memory); only on
        # the worker, so the task's command (its signature) stays the same
        cmd = ['env', 'ROOT_TTREECACHE_SIZE=0'] + cmd
    # results of the tasks the check may look at: all the tasks before it (the hadd of the -m2 run ranges compares with the
    # bdt_convert CV file, three levels up)
    need, todo = set(), list(t['deps'])
    while todo:
        d = todo.pop()
        if d not in need:
            need.add(d)
            todo.extend(by_name[d]['deps'])
    done_info = dict((n, {'info': state.data[n]['info']}) for n in need if state.data.get(n, {}).get('info') is not None)
    files['done.json'] = json.dumps(done_info)
    final = t['outputs'][0]
    dest = tmp_name(final) if t['intermediate'] else os.path.join(stage_dir, os.path.basename(final))
    sn = safe_name(t['name'])
    check = ' '.join(shlex.quote(x) for x in ['python3', '{code}/processing/run_production_det.py', '--worker-check', t['name'],
                                               '{out0}', 'done.json', '{out1}'] + worker_args(a))
    return dict(name=sn, cmd=['bash', '-c', ' '.join(shlex.quote(x) for x in cmd) + ' && ' + check],
                inputs=inputs, files=files, marker_dir=info_dir,
                outputs=[dict(name=os.path.basename(final), dest=dest), dict(name=sn + '.info.json', dest=os.path.join(info_dir, sn + '.info.json'))])


def copy_to(src, dst):
    '''copy src to dst through a temporary file, checking the size'''
    tmp = tmp_name(dst)
    remove(tmp)
    shutil.copyfile(src, tmp)
    if os.path.getsize(tmp) != os.path.getsize(src):
        remove(tmp)
        raise StepError('copy of %s to %s has the wrong size' % (src, dst))
    os.replace(tmp, dst)


def grid_runs(wiki, a, state, grids, g, run_dir):
    '''the run ranges of the CV files that have none yet, from grid jobs counting the events per run: True when all are known'''
    missing = [(k, p) for k, p in needed_cvs(wiki, a) if not grids.cached(k, p)]
    sub = []
    for key, path in missing:
        name = 'runs:cv_' + key
        st = state.data.get(name, {})
        out = os.path.join(run_dir, 'runs_cv_%s.json' % key)
        if st.get('status') == 'submitted':
            kind, marker = G.Grid.result(dict(marker_dir=run_dir, name=safe_name(name)), st.get('jobid'))
            if kind == 'done':
                counts = dict((int(r), n) for r, n in json.load(open(out)).items())
                grids.store(key, path, counts)
                st.update(status='done')
                state.save()
            elif kind == 'fail':
                st.update(status='failed', error='grid job failed, log %s' % os.path.join(run_dir, safe_name(name) + '.log'))
                state.save()
                log('%s: FAILED, log %s' % (name, os.path.join(run_dir, safe_name(name) + '.log')))
            continue
        if st.get('status') == 'failed' and st.get('tries', 0) >= a.max_tries:
            log('%s: failed %d times, not resubmitting (--reset it after a look at the log)' % (name, st['tries']))
            continue
        sub.append(dict(name=safe_name(name), marker_dir=run_dir,
                        cmd=['bash', '-c', 'python3 {code}/processing/run_production_det.py --worker-runs {in0} {out0}'],
                        inputs=[dict(path=path, mode='stream')], outputs=[dict(name=os.path.basename(out), dest=out)]))
    if sub:
        jobid = g.submit(sub, 'det_runs', memory=GRID_RES['runs'][0], lifetime=GRID_RES['runs'][1])
        if not a.dry_run:
            for t in sub:
                name = 'runs:' + t['name'][len('runs_'):]
                st = state.data.setdefault(name, {})
                st.update(status='submitted', jobid=jobid, tries=st.get('tries', 0) + 1)
            state.save()
    return not [k for k, p in needed_cvs(wiki, a) if not grids.cached(k, p)]


def run_grid(tasks, todo, by_name, consumers, a, state, g):
    '''one pass: collect finished jobs, submit the tasks that can run now; returns the number of jobs still running'''
    stage_dir = G.scratch_dir('det', 'final')
    info_dir = G.scratch_dir('det', 'jobs')
    queued = None
    for t in tasks:
        n = t['name']
        st = state.data.get(n, {})
        if st.get('status') != 'submitted':
            continue
        sn = safe_name(n)
        kind, marker = G.Grid.result(dict(marker_dir=info_dir, name=sn), st.get('jobid'))
        if kind == 'done':
            try:
                info = json.load(open(os.path.join(info_dir, sn + '.info.json')))
                final = t['outputs'][0]
                if t['intermediate']:
                    os.replace(tmp_name(final), final)
                else:
                    os.makedirs(os.path.dirname(final), exist_ok=True)
                    copy_to(os.path.join(stage_dir, os.path.basename(final)), final)
                    remove(os.path.join(stage_dir, os.path.basename(final)))
                with state.lock:
                    st.update(status='done', info=info, cmd=signature(t), seconds=round(marker.get('wall', 0)), error=None,
                              grid=dict((k, marker.get(k)) for k in ('host', 'wall', 'maxrss_mb', 'disk_mb')),
                              finished=time.strftime('%Y-%m-%d %H:%M:%S'), deleted=False)
                    state.save()
                todo.discard(n)
                log('%s: done on the grid (%.0f s, %.0f MB)' % (n, marker.get('wall', 0), marker.get('maxrss_mb', 0)))
            except Exception as e:
                with state.lock:
                    st.update(status='failed', error='after the grid job: %s' % e)
                    state.save()
                log('%s: FAILED after the grid job: %s' % (n, e))
        elif kind == 'fail':
            logf = os.path.join(info_dir, sn + '.log')
            with state.lock:
                st.update(status='failed', error='grid job failed at %s (rc %s), log %s' % (marker.get('stage'), marker.get('rc'), logf))
                state.save()
            log('%s: grid job FAILED at %s, log %s' % (n, marker.get('stage'), logf))
        else:
            if queued is None:
                queued = g.queued()
            if G.Grid.cluster(st.get('jobid')) not in queued:
                with state.lock:
                    st.update(status='failed', error='grid job %s left the queue without a result (jobsub_fetchlog for its log)' % st.get('jobid'))
                    state.save()
                log('%s: grid job %s lost' % (n, st.get('jobid')))
    delete_unneeded(tasks, a, state, todo, consumers)

    # tasks that can run now: to do, not submitted, every dependency done
    status = lambda n: state.data.get(n, {}).get('status')
    on_disk = detvar_on_disk(tasks, state, todo, consumers, lambda n: status(n) == 'submitted')
    ready = {}
    for t in tasks:
        n = t['name']
        if n not in todo or status(n) == 'submitted' or any(d in todo for d in t['deps']):
            continue
        if status(n) == 'failed' and state.data[n].get('tries', 0) >= a.max_tries:
            continue
        if starts_detvar(t):
            if on_disk >= a.max_detvar:
                continue
            on_disk += 1
        ready.setdefault(grid_kind(n), []).append(t)
    for kind, ts in ready.items():
        for i in range(0, len(ts), a.batch):
            batch = ts[i:i + a.batch]
            specs = [grid_task_of(t, a, by_name, state, stage_dir, info_dir) for t in batch]
            jobid = g.submit(specs, 'det_' + kind, memory=GRID_RES[kind][0], lifetime=GRID_RES[kind][1])
            if a.dry_run:
                log('%d %s tasks, e.g. %s: %s' % (len(batch), kind, batch[0]['name'], specs[0]['cmd'][2][:600]))
                continue
            with state.lock:
                for t in batch:
                    st = state.data.setdefault(t['name'], {})
                    st.update(status='submitted', jobid=jobid, tries=st.get('tries', 0) + 1, error=None, mode='grid',
                              submitted=time.strftime('%Y-%m-%d %H:%M:%S'))
                state.save()
    return sum(1 for t in tasks if status(t['name']) == 'submitted')


def write_det_cv_input(tasks, by_name, a, state):
    '''det_cv_input.txt of the -m2 files made so far (configurations/det_cv_input.txt format)'''
    lines = []
    for t in tasks:
        if t['name'].startswith('m3:') and state.data.get(t['name'], {}).get('status') == 'done':
            for d in t['deps'][1:]:
                lines.append('%s %s\n' % (by_name[d]['outputs'][0], t['outputs'][0]))
    with open(os.path.join(a.state_dir, 'det_cv_input.txt'), 'w') as f:
        f.write(''.join(lines) + 'end end\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-j', '--jobs', type=int, default=1, help='tasks run at the same time (default 1)')
    p.add_argument('--dry-run', action='store_true', help='print the tasks, run nothing')
    p.add_argument('--status', action='store_true', help='print the progress of each task and exit')
    p.add_argument('--reset', default=None, help='regular expression: forget the progress of the matching tasks, then exit')
    p.add_argument('--vars', default=None, help='comma separated: ' + ','.join(v for v, _ in VARS))
    p.add_argument('--runs', default=None, help='comma separated: ' + ','.join(RUNS))
    p.add_argument('--keep-work', action='store_true', help='do not delete the bdt_convert and run range files')
    p.add_argument('--chunk-events', type=int, default=CHUNK_EVENTS, help='CV events per run range (default %d)' % CHUNK_EVENTS)
    p.add_argument('--max-detvar', type=int, default=None,
                   help='DetVar files with pre-trim or bdt_convert files on disk at the same time (default: 10 with --grid, '
                   'otherwise -j, at least 2)')
    p.add_argument('--det-dir', default=os.path.join(OUT_DIR, 'DetVar'))
    p.add_argument('--work-dir', default=None, help='bdt_convert and run range files (default %s, with --grid %s)'
                   % (os.path.join(OUT_DIR, 'work', 'DetVar'), G.scratch_dir('det', 'work')))
    p.add_argument('--state-dir', default=os.path.join(HERE, 'state_det'), help='state file and logs (keep it off /pnfs)')
    p.add_argument('--bin-dir', default=os.path.join(TOP, 'bin'))
    p.add_argument('--wiki', default=os.path.join(TOP, 'numuCCana/mcc910wiki.txt'))
    p.add_argument('--config', default=os.path.join(HERE, 'config.txt'))
    p.add_argument('--trim-config', default=os.path.join(HERE, 'config_wc.txt'),
                   help='tree_trimmer config of the run ranges (config.txt with wcpselection T_pot as a POT tree)')
    p.add_argument('--pre-config', default=os.path.join(HERE, 'trimmer_config_det_pre.txt'),
                   help='tree_trimmer config of the pre-trim before bdt_convert (trees and branches kept)')
    p.add_argument('--weights', default=os.path.join(HERE, 'weights'), help='BDT weights directory of bdt_convert')
    p.add_argument('--chunks', default=None, help='run ranges file (default <state dir>/chunks.json)')
    p.add_argument('--grid', action='store_true', help='run the tasks as grid jobs (on a gpvm; no ROOT needed there)')
    p.add_argument('--wait', action='store_true', help='--grid: follow the jobs (and submit what becomes ready) until all are done')
    p.add_argument('--poll', type=int, default=300, help='--wait: seconds between checks (default 300)')
    p.add_argument('--max-tries', type=int, default=2, help='--grid: submissions of a task before giving up (default 2)')
    p.add_argument('--batch', type=int, default=500, help='--grid: tasks per submission (default 500)')
    p.add_argument('--app-dir', default=TOP, help='--grid: directory with bin/, lib/, lib64/ for the code tarball')
    p.add_argument('--worker-runs', nargs=2, metavar=('FILE', 'OUT'), help=argparse.SUPPRESS)
    p.add_argument('--worker-check', nargs=4, metavar=('TASK', 'FILE', 'DONE', 'OUT'), help=argparse.SUPPRESS)
    a = p.parse_args()
    a.state_dir = os.path.abspath(a.state_dir)
    a.run_dir = os.path.join(a.state_dir, 'run')
    a.chunks = a.chunks or os.path.join(a.state_dir, 'chunks.json')
    if a.work_dir is None:
        a.work_dir = G.scratch_dir('det', 'work') if a.grid else os.path.join(OUT_DIR, 'work', 'DetVar')
    a.max_detvar = a.max_detvar or (10 if a.grid else max(2, a.jobs))

    # grid worker: the run counts of an unprocessed CV file, or the check of a task's output (files read over xrootd)
    if a.worker_runs:
        os.environ['DET_GRID_WORKER'] = '1'
        counts = run_counts(a.worker_runs[0])
        with open(a.worker_runs[1], 'w') as f:
            json.dump(dict((str(r), n) for r, n in counts.items()), f)
        return 0
    if a.worker_check:
        os.environ['DET_GRID_WORKER'] = '1'
        name, output, done_file, info_file = a.worker_check
        tasks = build_tasks(rp.read_wiki(a.wiki), a, Grids(a.chunks, a.chunk_events, fixed=True).get)
        t = dict((x['name'], x) for x in tasks)[name]
        info = t['check'](output, json.load(open(done_file)))
        with open(info_file, 'w') as f:
            json.dump(info, f)
        log('%s: check passed' % name)
        return 0

    os.makedirs(os.path.join(a.state_dir, 'logs'), exist_ok=True)
    state = State(os.path.join(a.state_dir, 'state.json'), dry=a.dry_run or a.status)
    wiki = rp.read_wiki(a.wiki)
    grids = Grids(a.chunks, a.chunk_events)
    if not a.status and not a.dry_run:
        lockf = open(os.path.join(a.state_dir, '.lock'), 'w')
        try:
            fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            print('another run_production_det.py is using %s' % a.state_dir)
            return 1
    if a.grid and not a.status and not a.reset:
        # run ranges of the CV files: from grid jobs the first time (counting runs needs ROOT)
        g_runs = G.Grid(os.path.join(a.state_dir, 'grid'), a.app_dir, dry=a.dry_run, log=log, extra=[(a.wiki, 'processing/mcc910wiki.txt')])
        while not grid_runs(wiki, a, state, grids, g_runs, G.scratch_dir('det', 'runs')):
            if not a.wait or a.dry_run:
                log('the run ranges of the CV files come from grid jobs; run again when they are done (or use --wait)')
                return 0
            log('waiting for the run counts of the CV files, next check in %d s' % a.poll)
            time.sleep(a.poll)

    def get_grid(key, path):
        if not a.grid:
            return grids.get(key, path)
        c = grids.cached(key, path)
        if not c:
            raise StepError('run ranges of %s not known yet (run with --grid to get them)' % key)
        return c
    tasks = build_tasks(wiki, a, get_grid)
    by_name = dict((t['name'], t) for t in tasks)
    consumers = {}
    for t in tasks:
        for i in t['inputs']:
            consumers.setdefault(i, []).append(t)
    group = lambda n: re.sub(r'\.\d+$', '', n)  # task of a run range -> its task
    done = lambda t: state.data.get(t['name'], {}).get('status') == 'done' and state.data[t['name']].get('cmd') == signature(t)

    if a.status:
        ranges = {}
        for t in tasks:
            s = state.data.get(t['name'], {})
            st = s.get('status', 'pending') if (s.get('status') != 'done' or done(t)) else 'changed'
            err = ('  ' + s.get('error', '').strip().split('\n')[-1]) if st == 'failed' and s.get('error') else ''
            if group(t['name']) != t['name']:
                r = ranges.setdefault(group(t['name']), {'n': 0, 'err': ''})
                r['n'] += 1
                r[st] = r.get(st, 0) + 1
                r['err'] = r['err'] or err
                continue
            print('%-12s %-24s %s%s' % (st + ('*' if s.get('deleted') else ''), t['name'], os.path.basename(t['outputs'][0]), err))
        if ranges:
            print('\nrun ranges:')
            for g, r in ranges.items():
                print('%-24s %3d/%-3d done  %s%s' % (g, r.get('done', 0), r['n'], ', '.join('%d %s' % (r[k], k) for k in
                                                     ('running', 'failed', 'interrupted', 'changed') if r.get(k)), r['err']))
        print('(* output deleted after use)')
        return 0

    if a.reset:
        names = [n for n in list(state.data) if re.search(a.reset, n)]
        for n in names:
            state.data.pop(n)
        state.save()
        print('reset %d tasks: %s' % (len(names), ' '.join(names)))
        return 0

    # tasks to run: not done (or command / list changed), their dependents, and the producers of deleted inputs they need
    todo = set(t['name'] for t in tasks if not done(t) or (not t['intermediate'] and not all(os.path.exists(o) for o in t['outputs'])))
    changed = True
    while changed:
        changed = False
        for t in tasks:
            if t['name'] not in todo and any(d in todo and not done(by_name[d]) for d in t['deps']):
                todo.add(t['name'])
                changed = True
            if t['name'] in todo:
                for d in t['deps']:
                    if d not in todo and not all(os.path.exists(o) for o in by_name[d]['outputs']):
                        todo.add(d)
                        changed = True

    if a.grid:
        bad = [t['inputs'][0] for t in tasks if t['name'] in todo and t['name'].startswith('pre:') and not os.path.exists(t['inputs'][0])]
        if bad:
            print('inputs not found:\n  ' + '\n  '.join(bad))
            return 1
        g = G.Grid(os.path.join(a.state_dir, 'grid'), a.app_dir, weights=True, dry=a.dry_run, log=log,
                   extra=[(a.wiki, 'processing/mcc910wiki.txt'), (a.chunks, 'processing/det_chunks.json')])
        while True:
            n = run_grid(tasks, todo, by_name, consumers, a, state, g)
            if not a.wait or a.dry_run or not n:
                break
            log('%d grid jobs running, next check in %d s' % (n, a.poll))
            time.sleep(a.poll)
        if not a.dry_run:
            write_det_cv_input(tasks, by_name, a, state)
        left = [t['name'] for t in tasks if t['name'] in todo]
        log('%d of %d tasks done, %d grid jobs running%s' % (len(tasks) - len(left), len(tasks), n,
                                                              (', failed: ' + ' '.join(t for t in left if state.data.get(t, {}).get('status') == 'failed')) if left else ''))
        return 1 if left else 0

    if a.dry_run:
        shown = {}
        for t in tasks:
            g = group(t['name'])
            if g != t['name']:
                shown[g] = shown.get(g, 0) + 1
                if shown[g] > 1:
                    continue
            sig = signature(t)
            if t['cmd'][0] == 'hadd':
                sig = ' '.join(t['cmd'][:2] + [t['outputs'][0], '<the %d run range files>' % len(t['inputs'])])
            print('%-8s %-24s %s' % ('todo' if t['name'] in todo else 'done', t['name'], sig))
        print('\n(only the first run range of each task is shown: %s)' % ', '.join('%s %d' % (g, n) for g, n in shown.items()))
        print('%d tasks, %d to run' % (len(tasks), len(todo)))
        return 0

    host = socket.gethostname()
    for n, s in state.data.items():
        if s.get('status') == 'running' and s.get('pid'):
            if s.get('host') == host and pid_alive(s['pid']):
                print('%s is still running (pid %d) from an earlier run, wait for it or kill it' % (n, s['pid']))
                return 1
            if s.get('host') != host:
                print('WARNING: %s was running on %s (pid %d), make sure it is no longer running' % (n, s.get('host'), s['pid']))

    problems = []
    for x in ['bdt_convert', 'tree_trimmer']:
        if not os.access(os.path.join(a.bin_dir, x), os.X_OK):
            problems.append('no %s in %s' % (x, a.bin_dir))
    if subprocess.call('command -v hadd > /dev/null', shell=True) != 0:
        problems.append('no hadd (set up uboonecode)')
    if not os.access(os.path.join(a.bin_dir, 'merge_det'), os.X_OK):
        problems.append('no merge_det in %s' % a.bin_dir)
    else:
        out = subprocess.run([os.path.join(a.bin_dir, 'merge_det')], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.decode(errors='replace')
        if '-m3' not in out:
            problems.append('%s has no -m2 / -m3 (build it from the det_disk_saver branch)' % os.path.join(a.bin_dir, 'merge_det'))
    for f in [a.config, a.trim_config, a.pre_config, a.weights]:
        if not os.path.exists(f):
            problems.append('missing %s' % f)
    try:
        rp.root()
    except Exception as e:
        problems.append('no PyROOT: %s' % e)
    for t in tasks:
        if t['name'] in todo and t['name'].startswith('pre:') and not os.path.exists(t['inputs'][0]):
            problems.append('input of %s not found: %s' % (t['name'], t['inputs'][0]))
    if problems:
        print('\n'.join(problems))
        return 1

    for d in [a.det_dir, a.work_dir, a.run_dir]:
        os.makedirs(d, exist_ok=True)
    link = os.path.join(a.run_dir, 'weights')  # bdt_convert reads weights/*.xml from its working directory
    if os.path.realpath(link) != os.path.realpath(a.weights):
        remove(link)
        os.symlink(a.weights, link)

    def stop(signum, frame):
        if not rp.STOP.is_set():
            log('stopping (signal %d), killing the running apps' % signum)
        rp.STOP.set()
        for pid in list(rp.RUNNING):
            try:
                os.killpg(pid, signal.SIGTERM)
            except OSError:
                pass
    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)

    delete_unneeded(tasks, a, state, todo, consumers)  # files left by a run stopped before deleting them
    max_detvar = a.max_detvar

    def on_disk():
        return detvar_on_disk(tasks, state, todo, consumers, lambda n: n in running.values())

    failed = set()
    running = {}
    with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        while True:
            if not rp.STOP.is_set():
                for t in tasks:
                    if len(running) >= max(1, a.jobs):
                        break
                    n = t['name']
                    if n not in todo or n in running.values() or n in failed:
                        continue
                    if any(d in failed for d in t['deps']):
                        failed.add(n)
                        log('%s: skipped, needs %s' % (n, ', '.join(d for d in t['deps'] if d in failed)))
                        continue
                    if any(d in todo for d in t['deps']):
                        continue
                    if starts_detvar(t) and on_disk() >= max_detvar:
                        continue
                    running[pool.submit(run_task, t, a, state)] = n
            if not running:
                break
            fin, _ = wait(list(running), timeout=5, return_when=FIRST_COMPLETED)
            for fu in fin:
                n = running.pop(fu)
                if fu.result():
                    todo.discard(n)
                    delete_unneeded(tasks, a, state, todo, consumers)
                else:
                    failed.add(n)
            if rp.STOP.is_set() and not running:
                break

    write_det_cv_input(tasks, by_name, a, state)

    left = [t['name'] for t in tasks if t['name'] in todo]
    log('%d of %d tasks done%s; det_cv_input.txt in %s' % (len(tasks) - len(left), len(tasks),
                                                            (', not done: ' + ' '.join(left)) if left else '', a.state_dir))
    return 1 if left else 0


if __name__ == '__main__':
    sys.exit(main())
