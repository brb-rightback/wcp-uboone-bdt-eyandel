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
#   merge        ./bin/merge_hist -r0 -l0 -e2        -> its output in mcstat/0.log, the MC statistics TLee reads (run_stat.pl)
#                ./bin/merge_hist -r0 -l0 -e1        -> merge.root, the spectra TLee reads (only -e1 writes it; its output
#                                                     goes to the log, not to 0.log) (both after every conv)
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
# --breakdown: also the plot_hist breakdown plots in <work dir>/breakdown: the cv_input files whose cuts differ in the
# base folder's breakdown/file_ch.txt (the nu overlays, split by truth category) are converted again there with it
# (bconv:<n>), the others are linked from the work directory; then (breakdown) plot_hist <--breakdown-opts> -b1 -q1 there,
# by default -r0 -e3 -t2: the total uncertainty of the TLee covariance matrix (so after tlee) and the cross-section truth
# categories (-t2: XsecNumuCC0pinFV, XsecNumuCCNpinFV, XsecBkgCC, XsecNC, XsecCosmic of breakdown/file_ch.txt; the default
# -t0 stacks the LEE categories, which these files do not have) -> <channel>.png, canvas<n>.pdf,
# selection.pdf. On the grid its inputs are staged to scratch. plot_hist needs -q (batch, exits after the plots).
# --save: at the end, copy each configuration with a new finished tlee to numuCCana/hist_rootfiles_<config.txt>_<N>
# (--save-dir), in the layout of save_config.pl (also tlee.log, tlee_summary.txt, logs/); e.g. alone after a run:
#   python3 numuCCana/processing/run_validation.py --configs sp0p_N5_val --no-det --save
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
            'merge': ('2GB', '2h'), 'tlee': ('4GB', '8h'), 'plot': ('2GB', '2h')}
BREAKDOWN = 'breakdown'   # sub-directory of the work directory with the breakdown plots
# the tarballs of the grid jobs hold plain files only: extracted with the "data" filter where tarfile has it (no warning)
SAFE_TAR = {'filter': 'data'} if hasattr(tarfile, 'data_filter') else {}
# the last line read_TLee_v20 prints: after it, a non-zero exit code is ROOT's crash at exit (TTF::Cleanup, freetype),
# with every result already written
TLEE_FINISHED = ' ---> Finish all the program'


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


def breakdown_file_ch(base, config):
    '''the breakdown/file_ch.txt of a configuration (the deepest level that has one, as config_files), or None'''
    out, d = None, os.path.join(CONF, base)
    for level in [''] + config.split('/'):
        d = os.path.join(d, level) if level else d
        f = os.path.join(d, BREAKDOWN, 'file_ch.txt')
        if os.path.exists(f):
            out = open(f).read()
    return out


def file_ch_cuts(text):
    '''{input file: [its additional cuts]} of a file_ch.txt'''
    out = {}
    for l in text.split('\n'):
        t = l.split()
        if len(t) < 2 or t[0].startswith('#'):
            continue
        if t[0] == 'end':
            break
        out.setdefault(t[0], []).append(t[1])
    return out


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
    add(name='merge', kind='merge', cmd=['bash', '-c', './bin/merge_hist -r0 -l0 -e2 && ./bin/merge_hist -r0 -l0 -e1 1>&2'],
        out='merge.root', deps=list(convs), log_to='mcstat/0.log')
    add(name='tlee', kind='tlee', cmd=['./bin/read_TLee_v20'], out='file_collapsed_covariance_matrix.root',
        deps=['merge'] + stats + xfs + dets)
    if a.breakdown:
        # plot_hist breakdown plots in <work>/breakdown: the cv_input files whose file_ch cuts differ in breakdown/file_ch.txt
        # (the nu overlays, split by truth) converted again there with it, the others linked from the main conversion
        bfc = breakdown_file_ch(a.base, config)
        if bfc is None:
            raise rp.StepError('%s/%s: no %s/file_ch.txt' % (a.base, config, BREAKDOWN))
        cuts, bcuts = file_ch_cuts(files.get('file_ch.txt', '')), file_ch_cuts(bfc)
        deps = []
        for i, t in enumerate(table(files['cv_input.txt'])):
            if cuts.get(t[3]) == bcuts.get(t[3]):
                deps.append('conv:%02d' % i)
                continue
            add(name='bconv:%02d' % i, kind='conv', cmd=['./bin/convert_checkout_hist', t[3], t[4]], cwd=BREAKDOWN,
                out=os.path.join(BREAKDOWN, os.path.normpath(t[4])), deps=[], file=t[3], over={'file_ch.txt': bfc})
            deps.append('bconv:%02d' % i)
        opts = a.breakdown_opts.split()
        if '-e3' in opts:   # total uncertainty: the covariance matrix of TLee
            opts.append('-sfile_collapsed_covariance_matrix.root')
            deps.append('tlee')
        add(name='breakdown', kind='plot', cmd=['./bin/plot_hist'] + opts + ['-b1', '-q1'], cwd=BREAKDOWN,
            out=os.path.join(BREAKDOWN, 'selection.pdf'), deps=deps, over={'file_ch.txt': bfc})
    if a.no_det:
        files['TLee_config.txt'] = no_det(files['TLee_config.txt'])
    # signature: the command and the configuration files (TLee_config.txt only matters to tlee)
    h_all = hashlib.md5(''.join(k + files[k] for k in sorted(files) if k != 'TLee_config.txt').encode()).hexdigest()
    h_tlee = hashlib.md5(files['TLee_config.txt'].encode()).hexdigest()
    for t in tasks:
        # tlee, and the breakdown plots that draw its covariance matrix (-e3), are redone when TLee_config.txt changes (det)
        t['sig'] = ' '.join(t['cmd']) + ' | ' + h_all + (' ' + h_tlee if t['kind'] == 'tlee' or (t['kind'] == 'plot' and 'tlee' in t['deps']) else '')
        if t.get('over'):
            t['sig'] += ' ' + hashlib.md5(''.join(k + v for k, v in sorted(t['over'].items())).encode()).hexdigest()
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


def setup_breakdown(config, files, tasks, a):
    '''<work>/breakdown: the configuration with breakdown/file_ch.txt, DrawOption.cc (plot_hist runs it), and links to the
    inputs of plot_hist made in the work directory (the hist files not converted again, the TLee covariance matrix)'''
    plot = [t for t in tasks if t['kind'] == 'plot']
    if not plot:
        return
    w = work_dir(config, a)
    b = os.path.join(w, BREAKDOWN)
    for d in ('configurations', 'hist_rootfiles/XsFlux', 'hist_rootfiles/DetVar', 'mcstat'):
        os.makedirs(os.path.join(b, d), exist_ok=True)
    for name, text in dict(files, **plot[0]['over']).items():
        pth = os.path.join(b, 'configurations', name)
        if not os.path.exists(pth) or open(pth).read() != text:
            open(pth, 'w').write(text)
    shutil.copyfile(os.path.join(w, 'config.txt'), os.path.join(b, 'config.txt'))
    shutil.copyfile(os.path.join(rp.TOP, 'DrawOption.cc'), os.path.join(b, 'DrawOption.cc'))
    by = dict((t['name'], t) for t in tasks)
    links = [] if a.grid else [('bin', a.bin_dir)]
    for d in plot[0]['deps']:
        rel = by[d].get('final') or by[d]['out']
        if not rel.startswith(BREAKDOWN + '/'):
            links.append((rel, os.path.join(w, rel)))
    for rel, target in links:
        link = os.path.join(b, rel)
        target = os.path.relpath(target, os.path.dirname(link)) if rel != 'bin' else target
        if os.path.islink(link) and os.readlink(link) == target:
            continue
        remove(link)
        os.symlink(target, link)


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
        p = subprocess.Popen(t['cmd'], cwd=os.path.join(w, t.get('cwd', '')), stdout=out, stderr=lf if out is not lf else subprocess.STDOUT, start_new_session=True)
        rp.RUNNING[p.pid] = p
        rc = p.wait()
        rp.RUNNING.pop(p.pid, None)
        if out is not lf:
            out.close()
    o = os.path.join(w, t['out'])
    if rp.STOP.is_set():
        state.update(key, status='interrupted')
        return False
    if rc != 0 and t['kind'] == 'tlee' and TLEE_FINISHED in open(logf, errors='replace').read():
        log('%s: exit code %d after "%s" (ROOT crashes at exit), its results are complete' % (key, rc, TLEE_FINISHED.strip()))
        rc = 0
    if rc != 0 and t['kind'] == 'plot' and os.path.exists(o) and os.path.getsize(o) > 0:
        log('%s: exit code %d after writing %s (ROOT crashes at exit), kept' % (key, rc, t['out']))
        rc = 0
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
    for k, v in t.get('over', {}).items():   # breakdown/file_ch.txt, with the grid's URLs as the other files
        for pth in sorted(urls, key=len, reverse=True):
            v = v.replace(pth, urls[pth])
        files['configurations/' + k] = v
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
            if dt.get('log_to') and dt['log_to'] != rel:
                inputs.append(dict(path=scratch(config, dt['log_to']), mode='copy'))
                pre.append('ln -sf {in%d} %s' % (len(inputs) - 1, q(dt['log_to'])))
    if t['kind'] == 'plot':
        # its inputs are small and in the work directory (wherever they were made): staged to scratch, linked as in
        # <work>/breakdown; DrawOption.cc goes with the task
        files['DrawOption.cc'] = open(os.path.join(rp.TOP, 'DrawOption.cc')).read()
        w = work_dir(config, a)
        for d in t['deps']:
            rel = tasks_by_name[d].get('final') or tasks_by_name[d]['out']
            src, dst = os.path.join(w, rel), scratch(config, 'staged', rel)
            if not a.dry_run and (not os.path.exists(dst) or os.path.getsize(dst) != os.path.getsize(src)):
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                if os.path.exists(dst):
                    os.remove(dst)   # dCache files are not overwritten
                shutil.copyfile(src, dst)
            inputs.append(dict(path=dst, mode='copy'))
            pre.append('ln -sf {in%d} %s' % (len(inputs) - 1, q(os.path.relpath(rel, t['cwd']) if rel.startswith(t['cwd'] + '/') else rel)))
    cmd = ' '.join(q(urls.get(x, x)) for x in t['cmd'])   # the input files as in the worker's configuration
    outs = []
    if t['kind'] == 'tlee':   # everything it writes in the directory, as a tarball
        # ( ) not { }: grid_task.py fills the {in0}/{out0} placeholders with str.format
        run = ('touch .stamp && ( %s > .tlee_out 2>&1; rc=$?; cat .tlee_out; '
               'if [ $rc -ne 0 ]; then grep -q %s .tlee_out || exit $rc; echo "exit code $rc after the end (ROOT crash at exit), kept"; fi ) && '
               'find . -maxdepth 1 -type f -newer .stamp ! -name .stamp ! -name .tlee_out | tar czf {out0} -T -') % (cmd, q(TLEE_FINISHED))
        outs.append(dict(name='tlee_outputs.tar.gz', dest=scratch(config, 'tlee_outputs.tar.gz')))
    elif t['kind'] == 'plot':   # the plots, as a tarball; a crash at exit after selection.pdf is written is ignored
        run = ('touch .stamp && ( %s; rc=$?; [ -s selection.pdf ] || exit $rc ) && '
               'find . -maxdepth 1 -type f -newer .stamp ! -name .stamp | tar czf {out0} -T -') % cmd
        outs.append(dict(name='breakdown_plots.tar.gz', dest=scratch(config, 'breakdown_plots.tar.gz')))
    elif t.get('log_to') == t['out']:   # the output is what the command prints (merge: mcstat/0.log)
        run = '%s > %s && mv %s {out0}' % (cmd, q(t['out']), q(t['out']))
        outs.append(dict(name=os.path.basename(t['out']), dest=scratch(config, t['out'])))
    elif t.get('log_to'):
        run = '%s > %s && mv %s {out0} && mv %s {out1}' % (cmd, q(t['log_to']), q(t['out']), q(t['log_to']))
        outs += [dict(name=os.path.basename(t['out']), dest=scratch(config, t['out'])),
                 dict(name=os.path.basename(t['log_to']), dest=scratch(config, t['log_to']))]
    else:
        rel = t.get('final') or t['out']
        run = '%s && mv %s {out0}' % (cmd, q(os.path.relpath(t['out'], t.get('cwd') or '.')))
        outs.append(dict(name=os.path.basename(rel), dest=scratch(config, rel)))
    return dict(name=t['name'].replace(':', '_'), cmd=['bash', '-c', ' && '.join(pre) + ' && ' + run], inputs=inputs, files=files,
                outputs=outs, marker_dir=scratch(config, 'jobs'))


def collect(config, t, a):
    '''copy a finished grid task's outputs from scratch into the work directory'''
    w = work_dir(config, a)
    if t['kind'] == 'plot':
        with tarfile.open(scratch(config, 'breakdown_plots.tar.gz')) as tf:
            tf.extractall(os.path.join(w, t['cwd']), **SAFE_TAR)
        return
    if t['kind'] == 'tlee':
        tb = scratch(config, 'tlee_outputs.tar.gz')
        with tarfile.open(tb) as tf:
            tf.extractall(w, **SAFE_TAR)
        logf = os.path.join(scratch(config, 'jobs'), t['name'].replace(':', '_') + '.log')
        shutil.copyfile(logf, os.path.join(w, 'tlee.log'))
        return
    rels = [t.get('final') or t['out']] + ([t['log_to']] if t.get('log_to') and t['log_to'] != t['out'] else [])
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
            kind, marker = G.Grid.result(dict(marker_dir=scratch(config, 'jobs'), name=t['name'].replace(':', '_')), st.get('jobid'))
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


def save_breakdown(w, dest):
    '''the breakdown plots of the work directory w into dest/plots/breakdown (replacing those there), and the file_ch.txt
    they were made with as dest/configuration/breakdown_file_ch.txt'''
    bplots = glob.glob(os.path.join(w, BREAKDOWN, '*.png')) + glob.glob(os.path.join(w, BREAKDOWN, '*.pdf'))
    pd = os.path.join(dest, 'plots', BREAKDOWN)
    if os.path.exists(pd):
        shutil.rmtree(pd)
    os.makedirs(pd)
    for f in bplots:
        shutil.copy2(f, pd)
    os.makedirs(os.path.join(dest, 'configuration'), exist_ok=True)
    shutil.copy2(os.path.join(w, BREAKDOWN, 'configurations', 'file_ch.txt'), os.path.join(dest, 'configuration', 'breakdown_file_ch.txt'))
    return len(bplots)


def save(configs, a, state):
    '''--save: a copy of the work directory of every configuration with a finished tlee, in the layout of save_config.pl:
    <save dir>/hist_rootfiles_<config.txt>_<N> (N = the first free number) with the hist files, DetVar/, XsFlux/, mcstat/,
    merge.root, file_collapsed_covariance_matrix.root, plots/ (png, pdf, other root files), configuration/ (*.txt), and
    also tlee.log, tlee_summary.txt, logs/, and the breakdown plots (plots/breakdown/) if they are done. Copied (the work
    directory stays for reruns); each tlee result saved once; breakdown plots made after it are added to its folder.'''
    for config, (files, tasks, _, _) in configs.items():
        key = '%s:tlee' % config
        st = state.data.get(key, {})
        t = [x for x in tasks if x['name'] == 'tlee'][0]
        if not done(state, key, t):
            log('%s: tlee not done, not saved' % config)
            continue
        bkey = '%s:breakdown' % config
        bst = state.data.get(bkey, {})
        bt = [x for x in tasks if x['name'] == 'breakdown']
        bdone = done(state, bkey, bt[0]) if bt else bst.get('status') == 'done'
        w = work_dir(config, a)
        if st.get('saved_finished') == st.get('finished') and os.path.isdir(st.get('saved', '')):
            if bdone and st.get('saved_breakdown') != bst.get('finished'):
                if a.dry_run:
                    log('%s: would add the breakdown plots to %s' % (config, st['saved']))
                    continue
                n = save_breakdown(w, st['saved'])
                state.update(key, saved_breakdown=bst.get('finished'))
                log('%s: %d breakdown plots added to %s' % (config, n, st['saved']))
            else:
                log('%s: this tlee result is already saved in %s' % (config, st['saved']))
            continue
        name = open(os.path.join(w, 'config.txt')).read().strip()
        n = 0
        while os.path.exists(os.path.join(a.save_dir, 'hist_rootfiles_%s_%d' % (name, n))):
            n += 1
        dest = os.path.join(a.save_dir, 'hist_rootfiles_%s_%d' % (name, n))
        if a.dry_run:
            log('%s: would save to %s' % (config, dest))
            continue
        tmp = dest + '.tmp'
        if os.path.exists(tmp):
            shutil.rmtree(tmp)
        shutil.copytree(os.path.join(w, 'hist_rootfiles'), tmp)
        for d in ('mcstat', 'logs'):
            shutil.copytree(os.path.join(w, d), os.path.join(tmp, d))
        os.makedirs(os.path.join(tmp, 'plots'))
        os.makedirs(os.path.join(tmp, 'configuration'))
        for f in sorted(os.listdir(w)):
            p = os.path.join(w, f)
            if not os.path.isfile(p):
                continue
            if f in ('merge.root', 'file_collapsed_covariance_matrix.root', 'tlee.log', 'tlee_summary.txt', 'config.txt'):
                shutil.copy2(p, tmp)
            elif f.endswith(('.png', '.pdf', '.root')):
                shutil.copy2(p, os.path.join(tmp, 'plots'))
        for f in glob.glob(os.path.join(w, 'configurations', '*.txt')):
            shutil.copy2(f, os.path.join(tmp, 'configuration'))
        if bdone:
            save_breakdown(w, tmp)
        os.rename(tmp, dest)
        state.update(key, saved=dest, saved_finished=st.get('finished'), saved_breakdown=bst.get('finished') if bdone else None)
        log('%s: saved to %s' % (config, dest))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--configs', default='', help='comma separated configurations, folders of the base folder (e.g. sp0p_N5_val); '
                    'needed except with --reset')
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
    ap.add_argument('--save', action='store_true', help='at the end, copy the work directory of every configuration with a '
                    'new finished tlee to <save dir>/hist_rootfiles_<config.txt>_<N> (as save_config.pl)')
    ap.add_argument('--breakdown', action='store_true', help='also the plot_hist breakdown plots (truth categories of the nu '
                    'overlay, breakdown/file_ch.txt of the base folder) in <work dir>/breakdown')
    ap.add_argument('--breakdown-opts', default='-r0 -e3 -t2', help='--breakdown: plot_hist options (default "-r0 -e3 -t2": all '
                    'runs, total uncertainty from the TLee covariance matrix, the cross-section truth categories of breakdown/file_ch.txt '
                    '(numu CC 0p / Np in FV with the 45 MeV proton threshold, other CC, NC, cosmic, plus EXT and dirt); -b1 -q1 are added)')
    ap.add_argument('--save-dir', default=ANA, help='--save: where (default numuCCana, next to those of save_config.pl)')
    a = ap.parse_args()
    if not a.configs and not a.reset:
        ap.error('--configs is needed')
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
        print('forgot %d tasks: %s' % (len(names), ' '.join(names)))
        return 0

    for c, v in configs.items():
        setup(c, v[0], a)
        setup_breakdown(c, v[0], v[1], a)
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
        if a.save:
            save(configs, a, state)
        left = [(c, t['name']) for c, (_, tasks, _, _) in configs.items() for t in tasks if not done(state, '%s:%s' % (c, t['name']), t)]
        log('%d tasks not done%s' % (len(left), (': ' + ' '.join('%s:%s' % x for x in left[:20])) if left else ''))
        return 1 if left else 0

    if a.dry_run:
        for c, t in todo:
            print('%-22s %-10s (cd %s; %s%s)' % (c, t['name'], os.path.join(work_dir(c, a), t.get('cwd', '')), ' '.join(t['cmd']), ' > ' + t['log_to'] if t.get('log_to') else ''))
        print('%d tasks to run' % len(todo))
        if a.save:
            save(configs, a, state)
        return 0
    for x in ('convert_checkout_hist', 'merge_hist', 'read_TLee_v20', 'xf_cov_matrix', 'det_cov_matrix', 'stat_cov_matrix', 'stat_pred_cov_matrix') + (('plot_hist',) if a.breakdown else ()):
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
    if a.save:
        save(configs, a, state)
    log('%d tasks failed or not run%s' % (len(failed), (': ' + ' '.join(sorted(failed))) if failed else ''))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
