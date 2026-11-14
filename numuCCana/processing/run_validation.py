#!/usr/bin/env python3
# The model-validation tests of numuCCana configurations (cov_input.txt + TLee_config.txt), locally or on the grid, as the
# perl scripts of numuCCana do them by hand (setup_config.pl, convert_histo.pl, run_stat.pl, run_xf.pl, run_det.pl,
# read_TLee_v20), but each configuration in its own work directory, so that several run at once:
#
#   numuCCana/validation/<config with / as _>/
#       configurations/*.txt      the base files (e.g. configurations/opendata2/*.txt), overridden by those of every level
#                                 of the configuration's folder (as setup_config.pl), TLee_config.txt with --no-det: det 0
#       config.txt                <base>_<config> (run_stat.pl and the plots look at it)
#       bin -> <bin dir>          the apps (local mode)
#       hist_rootfiles/, hist_rootfiles/XsFlux/, hist_rootfiles/DetVar/, mcstat/, merge.root, the TLee plots
#       logs/<task>.log, tlee_summary.txt
#
# Tasks of a configuration (the perl scripts' commands, run from its work directory):
#   conv:<n>     ./bin/convert_checkout_hist <cv_input file> <its hist file>         one per cv_input.txt row
#   stat_data    ./bin/stat_cov_matrix -r0 -h1       -> mcstat/run_data_stat.root    configurations with "val" in the name
#   stat_pred    ./bin/stat_pred_cov_matrix -r0 -h1  -> mcstat/run_pred_stat.root     (as run_stat.pl)
#   xf:<N>       ./bin/xf_cov_matrix -r<N> -h1       -> hist_rootfiles/XsFlux/cov_<N>.root, N = 1 ... 17
#   (the xs_cov_matrix of run_xf.pl for "unfold" configurations is not done here: these are validation configurations)
#   det:<N>      ./bin/det_cov_matrix -r<N> -h1 -s7  -> the det_input.txt output of DetVar N (each N of det_input.txt)
#   merge        ./bin/merge_hist -r0 -l0 -e2        -> merge.root, its output in mcstat/0.log (after every conv)
#   tlee         ./bin/read_TLee_v20                 -> file_collapsed_covariance_matrix.root, plots; after all the others
# Then the "Summary of the TLee_config.txt tests" of read_TLee_v20 is read into tlee_summary.txt, and the table of every
# test of every configuration (with the description of the test from the comment above it in TLee_config.txt) is
# written to <validation dir>/summary.txt and summary.csv.
#
# Local (SL7 container with uboonecode set up):
#   python3 numuCCana/processing/run_validation.py --configs sp0p_N5_val --dry-run
#   python3 numuCCana/processing/run_validation.py --configs sp0p_N5_val,sp0p_N3_val -j 8 [--no-det]
# Grid (on a gpvm, natively; every task is a grid job, the files are read from /pnfs over xrootd: the /exp files of the
# configurations are replaced by their /pnfs copies, which must exist with the same size; merge and tlee copy their
# inputs to the worker; the outputs come back to the work directory):
#   python3 numuCCana/processing/run_validation.py --configs ... --grid --app-dir <build> [--wait]
# Both: --status, --reset <regex> (forget tasks, "<config>:<task>"), --summary (only rewrite the summary).
# --no-det: no detector systematics (no det tasks, det 0 in the work copy of TLee_config.txt); a later run without it
# redoes only tlee (and the det tasks). A task is redone when its command or the configuration files change.
import argparse, fcntl, glob, hashlib, json, os, re, shlex, shutil, signal, subprocess, sys, tarfile, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_production as rp
from run_production import State, log, remove
import grid as G

ANA = os.path.join(rp.TOP, 'numuCCana')
CONF = os.path.join(ANA, 'configurations')
# /exp areas with a /pnfs copy, for the grid
EXP_TO_PNFS = {'/exp/uboone/data/users/bbogart/numucc/': '/pnfs/uboone/persistent/users/bbogart/numucc/'}
XF_RUNS = list(range(1, 18))
DET_SEED = 7
GRID_RES = {'conv': ('2GB', '4h'), 'stat': ('2GB', '8h'), 'xf': ('3GB', '8h'), 'det': ('3GB', '8h'),
            'merge': ('2GB', '2h'), 'tlee': ('4GB', '8h')}


# ---------------------------------------------------------------------------------------------------------------------
# configurations
def config_files(base, config):
    '''{file name: text} of a configuration: the base folder's *.txt, overridden level by level (setup_config.pl)'''
    out = {}
    d = os.path.join(CONF, base)
    for level in [''] + config.split('/'):
        d = os.path.join(d, level) if level else d
        if not os.path.isdir(d):
            raise rp.StepError('no configuration folder %s' % d)
        for f in sorted(glob.glob(os.path.join(d, '*.txt'))):
            out[os.path.basename(f)] = open(f).read()
    for need in ('cov_input.txt', 'cv_input.txt', 'TLee_config.txt'):
        if need not in out:
            raise rp.StepError('%s/%s: no %s' % (base, config, need))
    return out


def table(text):
    '''rows of an input table (cv_input, xf_input, det_input): the lines before the first "-1" line, comments skipped'''
    rows = []
    for l in text.split('\n'):
        t = l.split()
        if not t or t[0].startswith('#'):
            continue
        if t[0] == '-1':
            break
        rows.append(t)
    return rows


def no_det(text):
    '''TLee_config.txt with det 0 in every test row (run index support target xsflux det add mcstat mcstatcor)'''
    out = []
    for l in text.split('\n'):
        t = l.split()
        if len(t) == 9 and not l.lstrip().startswith('#') and t[0].lstrip('-').isdigit():
            t[5] = '0'
            l = '  ' + ' '.join(t)
        out.append(l)
    return '\n'.join(out)


def to_grid(files, with_det=True):
    '''the configuration files with every input file replaced by the xrootd URL of its /pnfs copy (checked to exist with
    the same size; the detector-variation files only with_det); returns (files, problems, {path: URL})'''
    import grid_task
    paths = set()
    for name in ('cv_input.txt', 'xf_input.txt') + (('det_input.txt',) if with_det else ()):
        for t in table(files.get(name, '')):
            if len(t) > 3 and t[3].startswith('/'):
                paths.add(t[3])
    if with_det:
        for t in table(files.get('det_cv_input.txt', '').replace('end end', '-1')):
            paths.update(x for x in t if x.startswith('/'))
    sub, problems = {}, []
    for p in sorted(paths):
        q = p
        for e, pn in EXP_TO_PNFS.items():
            if p.startswith(e):
                q = pn + p[len(e):]
        if not q.startswith('/pnfs/'):
            problems.append('%s: no /pnfs copy known (EXP_TO_PNFS)' % p)
            continue
        if not os.path.exists(q):
            problems.append('%s: missing' % q)
            continue
        if q != p and os.path.exists(p) and os.path.getsize(p) != os.path.getsize(q):
            problems.append('%s and %s differ in size' % (p, q))
            continue
        sub[p] = grid_task.xrootd_url(q)
    out = {}
    for name, text in files.items():
        for p in sorted(sub, key=len, reverse=True):
            text = text.replace(p, sub[p])
        out[name] = text
    return out, problems, sub


def build_tasks(config, a):
    '''the tasks of a configuration: dict(name, kind, cmd (from the work directory), out (relative), deps, mcstat_from)'''
    files = config_files(a.base, config)
    tasks = []
    add = lambda **kw: tasks.append(kw)
    convs = []
    for i, t in enumerate(table(files['cv_input.txt'])):
        n = 'conv:%02d' % i
        add(name=n, kind='conv', cmd=['./bin/convert_checkout_hist', t[3], t[4]], out=os.path.normpath(t[4]), deps=[], file=t[3])
        convs.append(n)
    stats = []
    if 'val' in config:
        add(name='stat_data', kind='stat', cmd=['./bin/stat_cov_matrix', '-r0', '-h1'], out='hist_rootfiles/run_data_stat.root',
            final='mcstat/run_data_stat.root', deps=[])
        add(name='stat_pred', kind='stat', cmd=['./bin/stat_pred_cov_matrix', '-r0', '-h1'], out='hist_rootfiles/run_pred_stat.root',
            final='mcstat/run_pred_stat.root', deps=[])
        stats = ['stat_data', 'stat_pred']
    xfs = []
    if 'xf_input.txt' in files:
        for r in XF_RUNS:
            add(name='xf:%02d' % r, kind='xf', cmd=['./bin/xf_cov_matrix', '-r%d' % r, '-h1'], out='hist_rootfiles/XsFlux/cov_%d.root' % r, deps=[])
            xfs.append('xf:%02d' % r)
    dets = []
    if not a.no_det and 'det_input.txt' in files:
        outs = {}
        for t in table(files['det_input.txt']):
            outs.setdefault(int(t[2]), os.path.normpath(t[4]))
        for r, o in sorted(outs.items()):
            add(name='det:%02d' % r, kind='det', cmd=['./bin/det_cov_matrix', '-r%d' % r, '-h1', '-s%d' % DET_SEED], out=o, deps=[])
            dets.append('det:%02d' % r)
    add(name='merge', kind='merge', cmd=['./bin/merge_hist', '-r0', '-l0', '-e2'], out='merge.root', deps=list(convs), log_to='mcstat/0.log')
    add(name='tlee', kind='tlee', cmd=['./bin/read_TLee_v20'], out='file_collapsed_covariance_matrix.root',
        deps=['merge'] + stats + xfs + dets)
    if a.no_det:
        files['TLee_config.txt'] = no_det(files['TLee_config.txt'])
    # signature: the command and the configuration files (TLee_config.txt only matters to tlee)
    h_all = hashlib.md5(''.join(k + files[k] for k in sorted(files) if k != 'TLee_config.txt').encode()).hexdigest()
    h_tlee = hashlib.md5(files['TLee_config.txt'].encode()).hexdigest()
    for t in tasks:
        t['sig'] = ' '.join(t['cmd']) + ' | ' + h_all + (' ' + h_tlee if t['kind'] == 'tlee' else '')
    return files, tasks


def work_dir(config, a):
    return os.path.join(a.val_dir, config.replace('/', '_'))


def setup(config, files, a):
    '''the work directory of a configuration'''
    w = work_dir(config, a)
    for d in ('configurations', 'hist_rootfiles/XsFlux', 'hist_rootfiles/DetVar', 'mcstat', 'logs'):
        os.makedirs(os.path.join(w, d), exist_ok=True)
    for name, text in files.items():
        p = os.path.join(w, 'configurations', name)
        if not os.path.exists(p) or open(p).read() != text:
            open(p, 'w').write(text)
    open(os.path.join(w, 'config.txt'), 'w').write('%s_%s' % (a.base, config.replace('/', '_')))
    link = os.path.join(w, 'bin')
    if not a.grid and os.path.realpath(link) != os.path.realpath(a.bin_dir):
        remove(link)
        os.symlink(a.bin_dir, link)
    return w


# ---------------------------------------------------------------------------------------------------------------------
# local
def run_local(config, t, a, state):
    key = '%s:%s' % (config, t['name'])
    w = work_dir(config, a)
    logf = os.path.join(w, 'logs', t['name'].replace(':', '_') + '.log')
    t0 = time.time()
    state.update(key, status='running', sig=t['sig'], started=time.strftime('%Y-%m-%d %H:%M:%S'), error=None)
    log('%s: start' % key)
    for p in (t['out'], t.get('final')):
        if p:
            remove(os.path.join(w, p))
    with open(logf, 'w') as lf:
        lf.write(' '.join(t['cmd']) + '\n\n')
        lf.flush()
        out = open(os.path.join(w, t['log_to']), 'w') if t.get('log_to') else lf
        p = subprocess.Popen(t['cmd'], cwd=w, stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
        rp.RUNNING[p.pid] = p
        rc = p.wait()
        rp.RUNNING.pop(p.pid, None)
        if out is not lf:
            out.close()
    o = os.path.join(w, t['out'])
    if rp.STOP.is_set():
        state.update(key, status='interrupted')
        return False
    if rc != 0 or not os.path.exists(o) or os.path.getsize(o) == 0:
        state.update(key, status='failed', error='exit code %d%s, log %s' % (rc, '' if os.path.exists(o) else ', no %s' % t['out'], logf))
        log('%s: FAILED (exit code %d), see %s' % (key, rc, logf))
        return False
    if t.get('final'):
        os.replace(o, os.path.join(w, t['final']))
    if t['kind'] == 'tlee':
        shutil.copyfile(logf, os.path.join(w, 'tlee.log'))
    state.update(key, status='done', seconds=round(time.time() - t0), finished=time.strftime('%Y-%m-%d %H:%M:%S'))
    log('%s: done (%d s)' % (key, time.time() - t0))
    return True


# ---------------------------------------------------------------------------------------------------------------------
# grid
def scratch(config, *p):
    return G.scratch_dir('validation', config.replace('/', '_'), *p)


def grid_spec(config, t, gfiles, urls, tasks_by_name, a):
    '''the grid.py task: the configuration files written on the worker, the inputs of merge / tlee copied into place'''
    q = shlex.quote
    files = dict(('configurations/' + k, v) for k, v in gfiles.items())
    files['config.txt'] = open(os.path.join(work_dir(config, a), 'config.txt')).read()
    pre = ['mkdir -p hist_rootfiles/XsFlux hist_rootfiles/DetVar mcstat', 'ln -s {code}/bin bin']
    inputs = []
    if t['kind'] in ('merge', 'tlee'):
        # the outputs of the tasks before it, copied and linked to where the apps read them
        for d in t['deps']:
            dt = tasks_by_name[d]
            rel = dt.get('final') or dt['out']
            inputs.append(dict(path=scratch(config, rel), mode='copy'))
            pre.append('ln -sf {in%d} %s' % (len(inputs) - 1, q(rel)))
            if d == 'merge':
                inputs.append(dict(path=scratch(config, 'mcstat/0.log'), mode='copy'))
                pre.append('ln -sf {in%d} mcstat/0.log' % (len(inputs) - 1))
    cmd = ' '.join(q(urls.get(x, x)) for x in t['cmd'])   # the input files as in the worker's configuration
    outs = []
    if t['kind'] == 'tlee':   # everything it writes in the directory, as a tarball
        run = 'touch .stamp && %s && find . -maxdepth 1 -type f -newer .stamp ! -name .stamp | tar czf {out0} -T -' % cmd
        outs.append(dict(name='tlee_outputs.tar.gz', dest=scratch(config, 'tlee_outputs.tar.gz')))
    elif t.get('log_to'):
        run = '%s > %s && mv %s {out0} && mv %s {out1}' % (cmd, q(t['log_to']), q(t['out']), q(t['log_to']))
        outs += [dict(name=os.path.basename(t['out']), dest=scratch(config, t['out'])),
                 dict(name=os.path.basename(t['log_to']), dest=scratch(config, t['log_to']))]
    else:
        rel = t.get('final') or t['out']
        run = '%s && mv %s {out0}' % (cmd, q(t['out']))
        outs.append(dict(name=os.path.basename(rel), dest=scratch(config, rel)))
    return dict(name=t['name'].replace(':', '_'), cmd=['bash', '-c', ' && '.join(pre) + ' && ' + run], inputs=inputs, files=files,
                outputs=outs, marker_dir=scratch(config, 'jobs'))


def collect(config, t, a):
    '''copy a finished grid task's outputs from scratch into the work directory'''
    w = work_dir(config, a)
    if t['kind'] == 'tlee':
        tb = scratch(config, 'tlee_outputs.tar.gz')
        with tarfile.open(tb) as tf:
            tf.extractall(w)
        logf = os.path.join(scratch(config, 'jobs'), t['name'].replace(':', '_') + '.log')
        shutil.copyfile(logf, os.path.join(w, 'tlee.log'))
        return
    rels = [t.get('final') or t['out']] + ([t['log_to']] if t.get('log_to') else [])
    for rel in rels:
        dst = os.path.join(w, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(scratch(config, rel), dst + '.tmp')
        os.replace(dst + '.tmp', dst)


def run_grid(configs, a, state, g):
    '''one pass: collect finished jobs, submit the tasks that can run; returns the number still running'''
    queued = None
    ready = {}
    running = 0
    for config, (files, tasks, gfiles, _) in configs.items():
        by = dict((t['name'], t) for t in tasks)
        for t in tasks:
            key = '%s:%s' % (config, t['name'])
            st = state.data.get(key, {})
            if st.get('status') != 'submitted':
                continue
            kind, marker = G.Grid.result(dict(marker_dir=scratch(config, 'jobs'), name=t['name'].replace(':', '_')))
            if kind == 'done':
                try:
                    collect(config, t, a)
                    state.update(key, status='done', seconds=round(marker.get('wall', 0)), error=None, finished=time.strftime('%Y-%m-%d %H:%M:%S'),
                                 grid=dict((k, marker.get(k)) for k in ('host', 'wall', 'maxrss_mb', 'disk_mb')))
                    log('%s: done on the grid (%.0f s, %.0f MB)' % (key, marker.get('wall', 0), marker.get('maxrss_mb', 0)))
                except Exception as e:
                    state.update(key, status='failed', error='collecting the outputs: %s' % e)
                    log('%s: FAILED collecting the outputs: %s' % (key, e))
            elif kind == 'fail':
                logf = os.path.join(scratch(config, 'jobs'), t['name'].replace(':', '_') + '.log')
                state.update(key, status='failed', error='grid job failed at %s (rc %s), log %s' % (marker.get('stage'), marker.get('rc'), logf))
                log('%s: grid job FAILED, log %s' % (key, logf))
            else:
                if queued is None:
                    queued = g.queued()
                if G.Grid.cluster(st.get('jobid')) not in queued:
                    state.update(key, status='failed', error='grid job %s left the queue without a result' % st.get('jobid'))
                    log('%s: grid job lost' % key)
                else:
                    running += 1
        for t in tasks:
            key = '%s:%s' % (config, t['name'])
            if done(state, key, t) or state.data.get(key, {}).get('status') == 'submitted':
                continue
            if any(not done(state, '%s:%s' % (config, d), by[d]) for d in t['deps']):
                continue
            st = state.data.get(key, {})
            if st.get('status') == 'failed' and st.get('tries', 0) >= a.max_tries:
                continue
            ready.setdefault(t['kind'], []).append((config, t))
    for kind, items in ready.items():
        specs = [grid_spec(c, t, configs[c][2], configs[c][3], dict((x['name'], x) for x in configs[c][1]), a) for c, t in items]
        # one submission per configuration (task names repeat across configurations; markers are per configuration)
        for c in sorted(set(c for c, _ in items)):
            sub = [(s, t) for s, (cc, t) in zip(specs, items) if cc == c]
            jobid = g.submit([s for s, _ in sub], 'val_%s_%s' % (kind, c.replace('/', '_')), memory=GRID_RES[kind][0], lifetime=GRID_RES[kind][1])
            if a.dry_run:
                log('%s: %d %s tasks, e.g. %s' % (c, len(sub), kind, sub[0][0]['cmd'][2][:400]))
                continue
            for s, t in sub:
                key = '%s:%s' % (c, t['name'])
                st = state.data.get(key, {})
                state.update(key, status='submitted', jobid=jobid, sig=t['sig'], tries=st.get('tries', 0) + 1, error=None,
                             submitted=time.strftime('%Y-%m-%d %H:%M:%S'))
                running += 1
    return running


# ---------------------------------------------------------------------------------------------------------------------
def done(state, key, t):
    st = state.data.get(key, {})
    return st.get('status') == 'done' and st.get('sig') == t['sig']


def tlee_tests(text):
    '''{index: description} from TLee_config.txt (the comment line above each test row)'''
    out, desc = {}, ''
    for l in text.split('\n'):
        s = l.strip()
        if s.startswith('#'):
            desc = s.lstrip('#').strip()
            continue
        t = s.split()
        if len(t) == 9 and t[0].lstrip('-').isdigit():
            out[int(t[1])] = desc
    return out


def summarize(configs, a):
    '''the test table of read_TLee_v20 of every configuration -> <config>/tlee_summary.txt, <val dir>/summary.txt, .csv'''
    rows = []
    for config, (files, tasks, _, _) in configs.items():
        w = work_dir(config, a)
        logf = os.path.join(w, 'tlee.log')
        if not os.path.exists(logf):
            continue
        lines = open(logf, errors='replace').read().split('\n')
        try:
            i = [k for k, l in enumerate(lines) if 'Summary of the TLee_config.txt tests' in l][-1]
        except IndexError:
            continue
        desc = tlee_tests(files['TLee_config.txt'])
        own = []
        for l in lines[i + 2:]:
            t = l.split()
            if len(t) < 4 or not t[0].lstrip('-').isdigit():
                break
            idx, ndf, chi2, p = int(t[0]), int(t[1]), float(t[2]), float(t[3])
            r = dict(config=config, index=idx, ndf=ndf, chi2=chi2, p=p, test=desc.get(idx, ''))
            rows.append(r)
            own.append(r)
        open(os.path.join(w, 'tlee_summary.txt'), 'w').write(
            ''.join('%6d %5d %10.2f %10.3e  %s\n' % (r['index'], r['ndf'], r['chi2'], r['p'], r['test']) for r in own))
    if not rows:
        return
    hdr = '%-22s %6s %5s %10s %8s %10s  %s' % ('config', 'index', 'ndf', 'chi2', 'chi2/ndf', 'p-value', 'test')
    txt = [hdr, '-' * len(hdr)] + ['%-22s %6d %5d %10.2f %8.2f %10.3e  %s%s' % (
        r['config'], r['index'], r['ndf'], r['chi2'], r['chi2'] / r['ndf'] if r['ndf'] else 0, r['p'], r['test'], '   <-- p < 0.05' if r['p'] < 0.05 else '')
        for r in rows]
    open(os.path.join(a.val_dir, 'summary.txt'), 'w').write('\n'.join(txt) + '\n')
    with open(os.path.join(a.val_dir, 'summary.csv'), 'w') as f:
        f.write('config,index,ndf,chi2,p_value,test\n')
        for r in rows:
            f.write('%s,%d,%d,%.4f,%.6e,"%s"\n' % (r['config'], r['index'], r['ndf'], r['chi2'], r['p'], r['test'].replace('"', "'")))
    log('summary of %d tests in %s' % (len(rows), os.path.join(a.val_dir, 'summary.txt')))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--configs', required=True, help='comma separated configurations, folders of the base folder (e.g. sp0p_N5_val)')
    ap.add_argument('--base', default='opendata2', help='base folder in numuCCana/configurations (default opendata2)')
    ap.add_argument('--val-dir', default=os.path.join(ANA, 'validation'), help='work directories and summary (default numuCCana/validation)')
    ap.add_argument('--no-det', action='store_true', help='no detector systematics (det 0 in TLee_config.txt, no det tasks)')
    ap.add_argument('-j', '--jobs', type=int, default=4, help='local: tasks at the same time (default 4)')
    ap.add_argument('--bin-dir', default=os.path.join(rp.TOP, 'bin'), help='local: the apps')
    ap.add_argument('--grid', action='store_true', help='run the tasks as grid jobs (on a gpvm)')
    ap.add_argument('--app-dir', default=rp.TOP, help='--grid: directory with bin/, lib/, lib64/ for the code tarball')
    ap.add_argument('--wait', action='store_true', help='--grid: follow until everything is done')
    ap.add_argument('--poll', type=int, default=300)
    ap.add_argument('--max-tries', type=int, default=2)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--reset', default=None, help='regular expression on "<config>:<task>": forget those tasks, then exit')
    ap.add_argument('--summary', action='store_true', help='only rewrite the summary of the finished TLee runs')
    a = ap.parse_args()
    a.val_dir = os.path.abspath(a.val_dir)
    os.makedirs(a.val_dir, exist_ok=True)
    state = State(os.path.join(a.val_dir, 'state.json'), dry=a.dry_run or a.status)

    configs = {}
    for c in [x for x in a.configs.split(',') if x]:
        files, tasks = build_tasks(c, a)
        configs[c] = [files, tasks, None, None]   # + the grid files and their {path: URL}

    if a.status:
        for c, (files, tasks, _, _) in configs.items():
            for t in tasks:
                st = state.data.get('%s:%s' % (c, t['name']), {})
                s = st.get('status', 'pending') if (st.get('status') != 'done' or st.get('sig') == t['sig']) else 'changed'
                print('%-12s %-22s %-10s %s' % (s, c, t['name'], st.get('error') or ''))
        return 0
    if a.summary:
        summarize(configs, a)
        return 0
    lockf = open(os.path.join(a.val_dir, '.lock'), 'w')
    if not a.dry_run:
        try:
            fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            print('another run_validation.py is using %s' % a.val_dir)
            return 1
    if a.reset:
        names = [n for n in list(state.data) if re.search(a.reset, n)]
        for n in names:
            state.data.pop(n)
        state.save()
        print('forgot %d tasks' % len(names))
        return 0

    for c, v in configs.items():
        setup(c, v[0], a)
    todo = [(c, t) for c, (_, tasks, _, _) in configs.items() for t in tasks if not done(state, '%s:%s' % (c, t['name']), t)]

    if a.grid:
        problems = []
        for c, v in configs.items():
            v[2], p, v[3] = to_grid(v[0], with_det=not a.no_det)
            problems += ['%s: %s' % (c, x) for x in p]
        if problems:
            print('inputs not usable on the grid:\n  ' + '\n  '.join(problems))
            return 1
        g = G.Grid(os.path.join(a.val_dir, 'grid'), a.app_dir, dry=a.dry_run, log=log)
        while True:
            n = run_grid(configs, a, state, g)
            if not a.wait or a.dry_run or not n:
                break
            log('%d grid jobs running, next check in %d s' % (n, a.poll))
            time.sleep(a.poll)
        if not a.dry_run:
            summarize(configs, a)
        left = [(c, t['name']) for c, (_, tasks, _, _) in configs.items() for t in tasks if not done(state, '%s:%s' % (c, t['name']), t)]
        log('%d tasks not done%s' % (len(left), (': ' + ' '.join('%s:%s' % x for x in left[:20])) if left else ''))
        return 1 if left else 0

    if a.dry_run:
        for c, t in todo:
            print('%-22s %-10s (cd %s; %s%s)' % (c, t['name'], work_dir(c, a), ' '.join(t['cmd']), ' > ' + t['log_to'] if t.get('log_to') else ''))
        print('%d tasks to run' % len(todo))
        return 0
    for x in ('convert_checkout_hist', 'merge_hist', 'read_TLee_v20', 'xf_cov_matrix', 'det_cov_matrix', 'stat_cov_matrix', 'stat_pred_cov_matrix'):
        if not os.access(os.path.join(a.bin_dir, x), os.X_OK):
            print('no %s in %s' % (x, a.bin_dir))
            return 1

    def stop(signum, frame):
        log('stopping (signal %d), killing the running apps' % signum)
        rp.STOP.set()
        for pid in list(rp.RUNNING):
            try:
                os.killpg(pid, signal.SIGTERM)
            except OSError:
                pass
    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)
    failed, running = set(), {}
    from concurrent.futures import wait, FIRST_COMPLETED
    with ThreadPoolExecutor(max(1, a.jobs)) as pool:
        while True:
            for c, t in todo:
                key = '%s:%s' % (c, t['name'])
                if len(running) >= max(1, a.jobs) or rp.STOP.is_set():
                    break
                if key in running.values() or key in failed or done(state, key, t):
                    continue
                by = dict((x['name'], x) for x in configs[c][1])
                if any('%s:%s' % (c, d) in failed for d in t['deps']):
                    failed.add(key)
                    continue
                if not all(done(state, '%s:%s' % (c, d), by[d]) for d in t['deps']):
                    continue
                running[pool.submit(run_local, c, t, a, state)] = key
            if not running:
                break
            fin, _ = wait(list(running), timeout=5, return_when=FIRST_COMPLETED)
            for fu in fin:
                key = running.pop(fu)
                if not fu.result():
                    failed.add(key)
    summarize(configs, a)
    log('%d tasks failed or not run%s' % (len(failed), (': ' + ' '.join(sorted(failed))) if failed else ''))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
