#!/usr/bin/env python3
# merge_det in the original (-m1) format, CV and DetVar events in one file, for bdt_convert files too big to merge at
# once: the files are split in run ranges with tree_trimmer, merge_det runs on each range, the ranges are hadd'ed:
#
#   tree_trimmer <file> <range> <config> -o1 -k0 -r0 -y<lo> -z<hi>   the CV file and every DetVar file, run ranges of
#                                                                     about --chunk-events CV events
#   merge_det <cv range> <detvar range> <out range> -t<config>        per run range (the default mode, -m1; no -m option
#                                                                     is passed, so an older merge_det works too)
#   hadd <out> <out ranges>                                           per DetVar file: merge_det_<DetVar name>.root
#
# Exact: merge_det matches CV and DetVar events within a run and scales the POT per subrun, so the run ranges add up to
# the file merge_det would make from the whole files (the trees are in the order of the run ranges).
#
# Inputs: the CV bdt_convert file (the one with _nu_overlay_cv_ in its name, or --cv) and the DetVar bdt_convert files,
# default the run 4d BNB files below; or --files <list> (one path per line). Outputs in --out-dir (default: the directory of
# the CV file): merge_det_<DetVar file name without bdt_convert_><--suffix, default _v2>.root, so that the earlier
# merge_det files of the same name stay as they are; existing files not made by this script are never
# replaced. The run range files go to --work-dir and are deleted once used (--keep-work to keep them).
#
# Checks: every run range file has the events and POT subruns of its source in that range; every merge_det range has the
# same number of _cv and _det events and of T_pot_cv and T_pot_det subruns, no more events than either input range and,
# if both ranges have events, at least one match; every hadd file adds up to its ranges.
#
# Run inside the SL7 container with uboonecode set up (PyROOT, hadd):
#   python3 run_merge_det_ranges.py --dry-run
#   python3 run_merge_det_ranges.py -j 4
#   python3 run_merge_det_ranges.py --status
#   python3 run_merge_det_ranges.py --vars lya,sce -j 2          # a subset of the DetVar files (names as in the file)
#   python3 run_merge_det_ranges.py --reset 'merge:lya'           # forget the tasks matching a regex
# Resuming works as in run_production_det.py (whose run range grid and task runner are used here): every task writes to a
# temporary file renamed after its check, the progress is in <state dir>/state.json, finished tasks are skipped unless
# their command changed, and a deleted run range file is made again if a task that still has to run needs it.
import argparse, fcntl, os, re, signal, socket, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_production as rp
import run_production_det as pd
from run_production import StepError, State, log, remove, pid_alive

D4 = '/pnfs/uboone/persistent/users/uboonepro/surprise/detvar/BNB/run4d/'
DEFAULT_FILES = [D4 + 'bdt_convert_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_%s_surprise_reco2_hist_4d.root' % v
                 for v in ('cv', 'WMX', 'WMYZ', 'WMthetaXZ', 'WMthetaYZ', 'lya', 'lyd', 'lyr', 'recomb2', 'sce')]
CHUNK_EVENTS = 40000
# event trees checked when present (every file must have wcpselection/T_eval)
EVENT_TREES = ['wcpselection/T_eval', 'wcpselection/T_BDTvars', 'wcpselection/T_PFeval', 'wcpselection/T_KINEvars',
               'wcpselection/T_spacepoints', 'nuselection/NeutrinoSelectionFilter', 'lantern/EventTree']


def var_of(path):
    m = re.search(r'_nu_overlay_(.+?)_surprise_', os.path.basename(path))
    return m.group(1) if m else os.path.basename(path)[:-len('.root')]


def trees_in(path, suffix=''):
    '''the EVENT_TREES (with suffix) that exist in a file'''
    with rp.ROOT_LOCK:
        R = rp.root()
        f = R.TFile.Open(path)
        if not f or f.IsZombie():
            raise StepError('can not open %s' % path)
        try:
            return [t + suffix for t in EVENT_TREES if f.Get(t + suffix)]
        finally:
            f.Close()


def check_range(src, lo, hi):
    '''run range [lo, hi] of src: the events and POT subruns of src in that range (it can be empty)'''
    def check(path, done):
        trees = dict((t, []) for t in trees_in(path))
        if 'wcpselection/T_eval' not in trees:
            raise StepError('%s: no wcpselection/T_eval' % path)
        sel, psel = 'run>=%d && run<=%d' % (lo, hi), 'runNo>=%d && runNo<=%d' % (lo, hi)
        info = pd.query(path, dict(trees, **{'wcpselection/T_pot': ['runNo', 'subRunNo', 'pot_tor875']}),
                        counts={'n_out': ('wcpselection/T_eval', '!(%s)' % sel), 'n_pot_out': ('wcpselection/T_pot', '!(%s)' % psel)},
                        sums={'pot': ('wcpselection/T_pot', 'pot_tor875', '')})
        pd.same_entries(info, list(trees), path, empty=True)
        ref = pd.query(src, {}, counts={'n': ('wcpselection/T_eval', sel), 'n_pot': ('wcpselection/T_pot', psel)},
                       sums={'pot': ('wcpselection/T_pot', 'pot_tor875', psel)})
        if info['n_out'] or info['n_pot_out']:
            raise StepError('%s: %d events / %d subruns outside the runs %d-%d' % (path, info['n_out'], info['n_pot_out'], lo, hi))
        if info['wcpselection/T_eval'] != ref['n']:
            raise StepError('%s: %d events, %d in the runs %d-%d of %s' % (path, info['wcpselection/T_eval'], ref['n'], lo, hi, src))
        if info['wcpselection/T_pot'] != ref['n_pot'] or not pd.close(info['pot'], ref['pot']):
            raise StepError('%s: POT %g in %d subruns, %g in %d in the runs %d-%d of %s (is tree_trimmer built with the POT of '
                            '-y -z from the run limits?)' % (path, info['pot'], info['wcpselection/T_pot'], ref['pot'], ref['n_pot'], lo, hi, src))
        info['n'] = info['wcpselection/T_eval']
        return info
    return check


def check_merge(cv_task, det_task, parts=None):
    '''a merge_det (-m1) range (cv_task, det_task: its input ranges), or with parts the hadd of the ranges'''
    def check(path, done):
        cvt = trees_in(path, '_cv')
        det = trees_in(path, '_det')
        if 'wcpselection/T_eval_cv' not in cvt or 'wcpselection/T_eval_det' not in det:
            raise StepError('%s: no wcpselection/T_eval_cv / T_eval_det (not a merge_det -m1 file?)' % path)
        trees = dict((t, []) for t in cvt + det)
        info = pd.query(path, dict(trees, **{'wcpselection/T_pot_cv': ['runNo', 'subRunNo', 'pot_tor875'],
                                             'wcpselection/T_pot_det': ['runNo', 'subRunNo', 'pot_tor875']}),
                        sums={'pot_cv': ('wcpselection/T_pot_cv', 'pot_tor875', ''), 'pot_det': ('wcpselection/T_pot_det', 'pot_tor875', '')})
        # the Wire-Cell event trees of both sides have one entry per matched event (the trees of other directories are
        # only counted, and added up over the run ranges)
        pd.same_entries(info, [t for t in trees if t.startswith('wcpselection/')], path, empty=True)
        n, npot = info['wcpselection/T_eval_cv'], info['wcpselection/T_pot_cv']
        if npot != info['wcpselection/T_pot_det']:
            raise StepError('%s: %d T_pot_cv and %d T_pot_det subruns' % (path, npot, info['wcpselection/T_pot_det']))
        if parts:
            pd.check_sum(info, parts, list(trees) + ['wcpselection/T_pot_cv', 'wcpselection/T_pot_det', 'pot_cv', 'pot_det'], done, path)
            if n <= 0 or info['pot_cv'] <= 0:
                raise StepError('%s: %d matched events, CV POT %g' % (path, n, info['pot_cv']))
        else:
            ncv, ndet = done[cv_task]['info']['n'], done[det_task]['info']['n']
            if n > min(ncv, ndet):
                raise StepError('%s: %d matched events, the run ranges have %d CV and %d DetVar events' % (path, n, ncv, ndet))
            if ncv and ndet and n == 0:
                raise StepError('%s: no matched events from %d CV and %d DetVar events' % (path, ncv, ndet))
        return info
    return check


def build_tasks(a, cv, dets, grid):
    tasks, order = {}, []

    def add(name, **kw):
        kw.setdefault('content', None)
        kw.setdefault('intermediate', False)
        kw['name'] = name
        tasks[name] = kw
        order.append(name)

    chunk = lambda f, lo, hi: os.path.join(a.work_dir, 'ranges', '%s_%d_%d.root' % (os.path.basename(f)[:-len('.root')], lo, hi))
    trim = lambda src, lo, hi: [a.tree_trimmer, src, '{tmp}', a.trim_config, '-o1', '-k0', '-r0', '-y%d' % lo, '-z%d' % hi, '-v100000']
    rmin, rmax = grid['runs']
    for i, (lo, hi) in enumerate(grid['grid']):
        add('trim:cv.%02d' % i, cmd=trim(cv, lo, hi), inputs=[cv], outputs=[chunk(cv, lo, hi)], deps=[],
            check=check_range(cv, lo, hi), intermediate=True)
    for f in dets:
        v = var_of(f)
        parts = []
        for i, (lo, hi) in enumerate(grid['grid']):
            dlo, dhi = max(lo, rmin), min(hi, rmax)   # the DetVar ranges: those of the CV, within the CV runs
            ct, tt, mt = 'trim:cv.%02d' % i, 'trim:%s.%02d' % (v, i), 'merge:%s.%02d' % (v, i)
            add(tt, cmd=trim(f, dlo, dhi), inputs=[f], outputs=[chunk(f, dlo, dhi)], deps=[], check=check_range(f, dlo, dhi), intermediate=True)
            out = os.path.join(a.work_dir, 'ranges', 'merge_det_%s_%d_%d.root' % (os.path.basename(f)[len('bdt_convert_'):-len('.root')]
                                                                                  if os.path.basename(f).startswith('bdt_convert_') else
                                                                                  os.path.basename(f)[:-len('.root')], lo, hi))
            add(mt, cmd=[a.merge_det, tasks[ct]['outputs'][0], tasks[tt]['outputs'][0], '{tmp}', '-t' + a.config],
                inputs=[tasks[ct]['outputs'][0], tasks[tt]['outputs'][0]], outputs=[out], deps=[ct, tt], check=check_merge(ct, tt),
                intermediate=True)
            parts.append(mt)
        files = [tasks[t]['outputs'][0] for t in parts]
        add('merge:' + v, cmd=['hadd', '-f', '{tmp}'] + files, inputs=files, outputs=[final_name(a, f)], deps=parts,
            check=check_merge(None, None, parts=parts))
    return [tasks[n] for n in order]


def final_name(a, f):
    b = os.path.basename(f)
    b = b[len('bdt_convert_'):] if b.startswith('bdt_convert_') else b
    return os.path.join(a.out_dir, 'merge_det_' + b[:-len('.root')] + a.suffix + '.root')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-j', '--jobs', type=int, default=1, help='tasks run at the same time (default 1)')
    p.add_argument('--dry-run', action='store_true', help='print the tasks, run nothing')
    p.add_argument('--status', action='store_true', help='print the progress of each task and exit')
    p.add_argument('--reset', default=None, help='regular expression: forget the progress of the matching tasks, then exit')
    p.add_argument('--files', default=None, help='file with the bdt_convert files, one per line (default: the run 4d BNB DetVar files)')
    p.add_argument('--cv', default=None, help='the CV file (default: the file with _nu_overlay_cv_ in its name)')
    p.add_argument('--vars', default=None, help='comma separated DetVar names (as in the file names, e.g. lya,WMX): only these')
    p.add_argument('--keep-work', action='store_true', help='do not delete the run range files')
    p.add_argument('--chunk-events', type=int, default=CHUNK_EVENTS, help='CV events per run range (default %d)' % CHUNK_EVENTS)
    p.add_argument('--out-dir', default=None, help='output directory (default: the directory of the CV file)')
    p.add_argument('--suffix', default='_v2', help='added to the output names before .root (default _v2; "" for none)')
    p.add_argument('--work-dir', default='/pnfs/uboone/scratch/users/%s/merge_det_ranges_work' % os.environ.get('USER', 'user'))
    p.add_argument('--state-dir', default=os.path.join(os.getcwd(), 'state_merge_det_ranges'), help='state file and logs (keep it off /pnfs)')
    p.add_argument('--bin-dir', default=os.path.join(rp.TOP, 'bin'))
    p.add_argument('--merge-det', default=None, help='the merge_det to use (default <bin dir>/merge_det)')
    p.add_argument('--tree-trimmer', default=None, help='the tree_trimmer to use (default <bin dir>/tree_trimmer)')
    p.add_argument('--config', default=os.path.join(HERE, 'config.txt'), help='merge_det -t config')
    p.add_argument('--trim-config', default=os.path.join(HERE, 'config_wc.txt'),
                   help='tree_trimmer config of the run ranges (config.txt with wcpselection T_pot as a POT tree)')
    a = p.parse_args()
    a.merge_det = a.merge_det or os.path.join(a.bin_dir, 'merge_det')
    a.tree_trimmer = a.tree_trimmer or os.path.join(a.bin_dir, 'tree_trimmer')
    a.state_dir = os.path.abspath(a.state_dir)
    a.run_dir = os.path.join(a.state_dir, 'run')
    os.makedirs(os.path.join(a.state_dir, 'logs'), exist_ok=True)
    state = State(os.path.join(a.state_dir, 'state.json'), dry=a.dry_run or a.status)

    files = [l.strip() for l in open(a.files) if l.strip() and not l.startswith('#')] if a.files else DEFAULT_FILES
    cv = a.cv or next((f for f in files if '_nu_overlay_cv_' in os.path.basename(f)), None)
    if cv is None:
        print('no CV file (give --cv)')
        return 1
    dets = [f for f in files if f != cv and (not a.vars or var_of(f) in a.vars.split(','))]
    if len(set(var_of(f) for f in dets)) != len(dets):
        print('two DetVar files with the same name: %s' % ' '.join(var_of(f) for f in dets))
        return 1
    a.out_dir = a.out_dir or os.path.dirname(cv)
    missing = [f for f in [cv] + dets if not os.path.exists(f)]
    if missing and not a.status:
        print('not found:\n  ' + '\n  '.join(missing))
        return 1
    grids = pd.Grids(os.path.join(a.state_dir, 'chunks.json'), a.chunk_events)
    tasks = build_tasks(a, cv, dets, grids.get('cv', cv) if not missing else dict(grid=[[0, pd.RUN_MAX]], runs=[0, pd.RUN_MAX]))
    by_name = dict((t['name'], t) for t in tasks)
    consumers = {}
    for t in tasks:
        for i in t['inputs']:
            consumers.setdefault(i, []).append(t)
    group = lambda n: re.sub(r'\.\d+$', '', n)
    done = lambda t: state.data.get(t['name'], {}).get('status') == 'done' and state.data[t['name']].get('cmd') == pd.signature(t)

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
            print('%-12s %-22s %s%s' % (st, t['name'], t['outputs'][0], err))
        print('\nrun ranges:')
        for g, r in ranges.items():
            print('%-22s %3d/%-3d done  %s%s' % (g, r.get('done', 0), r['n'], ', '.join('%d %s' % (r[k], k) for k in
                                                 ('running', 'failed', 'interrupted', 'changed') if r.get(k)), r['err']))
        return 0

    lockf = open(os.path.join(a.state_dir, '.lock'), 'w')
    if not a.dry_run:
        try:
            fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            print('another run_merge_det_ranges.py is using %s' % a.state_dir)
            return 1
    if a.reset:
        names = [n for n in list(state.data) if re.search(a.reset, n)]
        for n in names:
            state.data.pop(n)
        state.save()
        print('reset %d tasks: %s' % (len(names), ' '.join(names)))
        return 0

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

    # never replace a final file this script did not make
    foreign = [t['outputs'][0] for t in tasks if not t['intermediate'] and t['name'] in todo and os.path.exists(t['outputs'][0])
               and state.data.get(t['name'], {}).get('status') != 'done']
    if a.dry_run:
        shown = {}
        for t in tasks:
            g = group(t['name'])
            if g != t['name']:
                shown[g] = shown.get(g, 0) + 1
                if shown[g] > 1:
                    continue
            sig = pd.signature(t)
            if t['cmd'][0] == 'hadd':
                sig = ' '.join(t['cmd'][:2] + [t['outputs'][0], '<the %d run range files>' % len(t['inputs'])])
            print('%-5s %-22s %s' % ('todo' if t['name'] in todo else 'done', t['name'], sig))
        print('\nCV %s: %d run ranges of about %d CV events' % (cv, sum(1 for t in tasks if t['name'].startswith('trim:cv.')), a.chunk_events))
        print('(only the first run range of each task is shown: %s)' % ', '.join('%s %d' % (g, n) for g, n in shown.items()))
        if foreign:
            print('\nEXIST and were not made by this script (the run would stop; remove them or use --out-dir):\n  ' + '\n  '.join(foreign))
        print('%d tasks, %d to run' % (len(tasks), len(todo)))
        return 0
    if foreign:
        print('these outputs exist and were not made by this script, not replacing them (remove them or use --out-dir):\n  ' + '\n  '.join(foreign))
        return 1

    host = socket.gethostname()
    for n, s in state.data.items():
        if s.get('status') == 'running' and s.get('pid'):
            if s.get('host') == host and pid_alive(s['pid']):
                print('%s is still running (pid %d) from an earlier run, wait for it or kill it' % (n, s['pid']))
                return 1
            if s.get('host') != host:
                print('WARNING: %s was running on %s (pid %d), make sure it is no longer running' % (n, s.get('host'), s['pid']))
    problems = []
    for x in (a.merge_det, a.tree_trimmer):
        if not os.access(x, os.X_OK):
            problems.append('not executable: %s' % x)
    if subprocess.call('command -v hadd > /dev/null', shell=True) != 0:
        problems.append('no hadd (set up uboonecode)')
    for f in (a.config, a.trim_config):
        if not os.path.exists(f):
            problems.append('missing %s' % f)
    try:
        rp.root()
    except Exception as e:
        problems.append('no PyROOT: %s' % e)
    if problems:
        print('\n'.join(problems))
        return 1
    for d in (a.out_dir, a.work_dir, a.run_dir):
        os.makedirs(d, exist_ok=True)

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

    pd.delete_unneeded(tasks, a, state, todo, consumers)
    failed, running = set(), {}
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
                    running[pool.submit(pd.run_task, t, a, state)] = n
            if not running:
                break
            fin, _ = wait(list(running), timeout=5, return_when=FIRST_COMPLETED)
            for fu in fin:
                n = running.pop(fu)
                if fu.result():
                    todo.discard(n)
                    pd.delete_unneeded(tasks, a, state, todo, consumers)
                else:
                    failed.add(n)
            if rp.STOP.is_set() and not running:
                break
    left = [t['name'] for t in tasks if t['name'] in todo]
    log('%d of %d tasks done%s' % (len(tasks) - len(left), len(tasks), (', not done: ' + ' '.join(left)) if left else ''))
    return 1 if left else 0


if __name__ == '__main__':
    sys.exit(main())
