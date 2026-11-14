#!/usr/bin/env python3
# The numuCC checkouts of the BNB beam-on data from the POT-trimmed _v2 files of run_trim_v2.py: tree_trimmer with the
# numuCC tree config, the flags of run_production.py's final trim of data (-o1 -k1 -c1: duplicates removed, Lantern
# failures removed, default -s), for the 13 files in INPUTS (runs 1-3, the run 1 open trigger period, the retuple runs
# 4a-5 and the open data of runs 1, 3, 4a, 4b):
#
#   <uboonepro surprise area>/.../trim_checkout_<name>_v2.root  ->  checkout_<name>_v2.root
#                                                                  in /pnfs/uboone/persistent/users/bbogart/numucc/
#                                                                  and /exp/uboone/data/users/bbogart/numucc/
#
# Per file: trim, then file_info.py of input and output (tree entries, T_pot POT), then the check, then the copies.
# Check: every event tree has the same entries, and the events, T_pot and SubRun entries and POT are exactly those of the
# input (nothing is expected to go: no run cuts with the default -s, no Lantern failures left in the _v2 files, no
# duplicates expected); --max-loss allows a fraction to go after a look. The input files are never changed.
# Each output is written to a temporary name and renamed after its check; nothing that exists is replaced unless this
# script made it.
#
# Local (SL7 container with uboonecode set up: tree_trimmer and PyROOT):
#   python3 numuCCana/processing/run_trim_numucc.py --dry-run
#   python3 numuCCana/processing/run_trim_numucc.py -j 4 [--bin-dir <build>/bin]
# Grid (on a gpvm, no ROOT needed; the apps come from the code tarball of --app-dir, which must have the TFile::Open
# version of tree_trimmer, as the input is streamed over xrootd):
#   python3 numuCCana/processing/run_trim_numucc.py --grid --app-dir <build>       # submit what is not done
#   python3 numuCCana/processing/run_trim_numucc.py --grid --app-dir <build> --wait # submit, then follow until done
#   python3 numuCCana/processing/run_trim_numucc.py --grid                         # follow: check and copy finished jobs
#   (a failed or lost job is resubmitted on the next --grid run, up to --max-tries submissions per file)
# Both:
#   python3 numuCCana/processing/run_trim_numucc.py --status
#   python3 numuCCana/processing/run_trim_numucc.py --reset <regex>    # forget the progress of matching files
#   python3 numuCCana/processing/run_trim_numucc.py --delete-old       # remove the old checkout_<name>.root (data,
#                                                                      # /pnfs and /exp) of files that are done; asks
# Progress in <state dir>/state.json, logs in <state dir>/logs/.
import argparse, fcntl, json, os, re, shutil, signal, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_production as rp
from run_production import State, log, tmp_name, remove
import grid as G

V2_DIR = '/pnfs/uboone/persistent/users/uboonepro/surprise/'
INPUTS = [
    'run1_full_samples/wc_processed/BNB/trim_checkout_MCC9.10_Run123_v10_04_07_20_BNB_beam_on_data_surprise_reco2_hist_1_v2.root',
    'run2_full_samples/wc_processed/BNB/trim_checkout_MCC9.10_Run123_v10_04_07_20_BNB_beam_on_data_surprise_reco2_hist_2_v2.root',
    'run3_full_samples/wc_processed/BNB/trim_checkout_MCC9.10_Run123_v10_04_07_20_BNB_beam_on_data_surprise_reco2_hist_3_v2.root',
    'run1_full_samples/wc_processed/BNB/trim_checkout_MCC9.10_Run1_OpenTrig_v10_04_07_26_BNB_beam_on_data_surprise_reco2_hist_v2.root',
    'opendata/BNB/trim_checkout_MCC9.10_Run123_v10_04_07_20_BNB_beam_on_data_surprise_reco2_hist_1_5e19opendata_v2.root',
    'opendata/BNB/trim_checkout_MCC9.10_Run123_v10_04_07_20_BNB_beam_on_data_surprise_reco2_hist_3_1e19opendata_v2.root',
    'retuple/BNB/trim_checkout_MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_4a_v2.root',
    'retuple/BNB/trim_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_beam_on_metapatch_retuple_retuple_hist_v2.root',
    'retuple/BNB/trim_checkout_MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_4c_v2.root',
    'retuple/BNB/trim_checkout_MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_4d_v2.root',
    'retuple/BNB/trim_checkout_MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_5_v2.root',
    'retuple/BNB/trim_checkout_MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_opendata_19550_v2.root',
    'retuple/BNB/trim_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_beam_on_metapatch_retuple_retuple_hist_opendata_20700_v2.root',
]
PNFS_OUT = '/pnfs/uboone/persistent/users/bbogart/numucc'
EXP_OUT = '/exp/uboone/data/users/bbogart/numucc'
FLAGS = ['-o1', '-k1', '-c1']
GRID_RES = dict(memory='2GB', lifetime='4h')


def out_name(inp):
    '''trim_checkout_<name>_v2.root -> checkout_<name>_v2.root'''
    b = os.path.basename(inp)
    if not (b.startswith('trim_checkout_') and b.endswith('_v2.root')):
        raise ValueError(inp)
    return b[len('trim_'):]


def old_name(inp):
    '''the checkout of the same data made by run_production.py from the unprocessed file (none for the open trigger)'''
    return out_name(inp)[:-len('_v2.root')] + '.root'


def samples():
    return [dict(input=os.path.join(V2_DIR, i), name=out_name(i), pnfs=os.path.join(PNFS_OUT, out_name(i)),
                 exp=os.path.join(EXP_OUT, out_name(i)), old=[os.path.join(d, old_name(i)) for d in (PNFS_OUT, EXP_OUT)])
            for i in INPUTS]


# ---------------------------------------------------------------------------------------------------------------------
def check(s, info_in, info_out, max_loss):
    '''the checks of the trimmed file against its input; returns a summary. Nothing should be lost (no run cuts with the
    default -s, the _v2 files have no Lantern failures left, duplicates are not expected): events, T_pot and SubRun
    entries and POT must be those of the input, unless --max-loss allows a fraction of each to go'''
    rp.check_events(info_out, s['name'])
    ev = 'wcpselection/T_eval'
    for key, what in ((ev, 'events'), (rp.POT_TREE[0], 'T_pot entries'), ('nuselection/SubRun', 'SubRun entries'), ('pot', 'T_pot POT')):
        x_in, x_out = info_in[key], info_out[key]
        tol = 1e-6 * x_in if key == 'pot' else 0
        if x_out > x_in + tol:
            raise rp.StepError('%s: %s %s above the input %s' % (s['name'], what, x_out, x_in))
        if x_out < x_in * (1 - max_loss) - tol:
            raise rp.StepError('%s: %s %s, the input has %s (%.3g%% lost; --max-loss %g)'
                               % (s['name'], what, x_out, x_in, 100. * (x_in - x_out) / x_in if x_in else 0, max_loss))
    return dict(events_in=info_in[ev], events_out=info_out[ev], subruns_in=info_in['nuselection/SubRun'], subruns_out=info_out['nuselection/SubRun'],
                tpot_in=info_in[rp.POT_TREE[0]], tpot_out=info_out[rp.POT_TREE[0]], pot_in=info_in['pot'], pot_out=info_out['pot'])


def copy_to(src, dst):
    '''copy src to dst through a temporary file, checking the size'''
    tmp = tmp_name(dst)
    remove(tmp)
    shutil.copyfile(src, tmp)
    if os.path.getsize(tmp) != os.path.getsize(src):
        remove(tmp)
        raise rp.StepError('copy of %s to %s has the wrong size' % (src, dst))
    os.replace(tmp, dst)


def place(s, produced, state, keep_produced):
    '''produced (a checked file: the local temporary or the grid scratch output) -> the /pnfs and /exp outputs'''
    if os.path.abspath(produced) == os.path.abspath(tmp_name(s['pnfs'])):
        os.replace(produced, s['pnfs'])
    else:
        copy_to(produced, s['pnfs'])
        if not keep_produced:
            remove(produced)
    copy_to(s['pnfs'], s['exp'])
    state.update(s['name'], status='done', size=os.path.getsize(s['pnfs']), finished=time.strftime('%Y-%m-%d %H:%M:%S'))
    log('%s: done (%s, %s)' % (s['name'], s['pnfs'], s['exp']))


def finish(s, produced, info_file, state, a):
    '''check a produced file from its file_info.py json and move it on; returns True if done'''
    try:
        (p_in, info_in), (p_out, info_out) = json.load(open(info_file))
        summary = check(s, info_in, info_out, a.max_loss)
        state.update(s['name'], check=summary)
        log('%s: %d -> %d events, %d -> %d subruns, T_pot POT %.4g -> %.4g' % (
            s['name'], summary['events_in'], summary['events_out'], summary['subruns_in'], summary['subruns_out'],
            summary['pot_in'], summary['pot_out']))
        place(s, produced, state, keep_produced=False)
        return True
    except Exception as e:   # a check failure needs a look: not resubmitted automatically (--reset to retry)
        state.update(s['name'], status='failed', error=str(e), no_retry=True)
        log('%s: FAILED %s (not resubmitted; --reset it to try again)' % (s['name'], e))
        return False


# ---------------------------------------------------------------------------------------------------------------------
# local
def run_local(s, a, state):
    tmp = tmp_name(s['pnfs'])
    info_file = os.path.join(a.state_dir, 'info', s['name'] + '.json')
    logfile = os.path.join(a.state_dir, 'logs', s['name'] + '.log')
    cmds = [[os.path.join(a.bin_dir, 'tree_trimmer'), s['input'], tmp, a.config] + FLAGS,
            [sys.executable, os.path.join(HERE, 'file_info.py'), s['input'], tmp]]
    state.update(s['name'], status='running', mode='local', cmd=' '.join(cmds[0]), started=time.strftime('%Y-%m-%d %H:%M:%S'),
                 log=logfile, error=None)
    log('%s: start' % s['name'])
    remove(tmp)
    with open(logfile, 'w') as lf:
        for k, cmd in enumerate(cmds):
            if rp.STOP.is_set():
                return False
            lf.write(' '.join(cmd) + '\n\n')
            lf.flush()
            out = open(info_file, 'w') if k == 1 else lf
            p = subprocess.Popen(cmd, stdout=out, stderr=lf, start_new_session=True)
            rp.RUNNING[p.pid] = p
            rc = p.wait()
            rp.RUNNING.pop(p.pid, None)
            if out is not lf:
                out.close()
            if rc != 0:
                state.update(s['name'], status='failed', error='%s exit code %d (log %s)' % (os.path.basename(cmd[0]), rc, logfile))
                log('%s: FAILED, see %s' % (s['name'], logfile))
                remove(tmp)
                return False
    return finish(s, tmp, info_file, state, a)


# ---------------------------------------------------------------------------------------------------------------------
# grid
def grid_paths(s):
    d = G.scratch_dir('trim_numucc')
    return os.path.join(d, s['name']), os.path.join(d, s['name'][:-len('.root')] + '_info.json')


def grid_task(s):
    out, info = grid_paths(s)
    return dict(name=s['name'][:-len('.root')],
                cmd=['bash', '-c', '{code}/bin/tree_trimmer {in0} {out0} {code}/processing/%s %s && '
                     'python3 {code}/processing/file_info.py {in0} {out0} > {out1}' % (os.path.basename(GRID_CONFIG), ' '.join(FLAGS))],
                inputs=[dict(path=s['input'], mode='stream')],
                outputs=[dict(name=s['name'], dest=out), dict(name=os.path.basename(info), dest=info)],
                marker_dir=os.path.dirname(out))


GRID_CONFIG = os.path.join(HERE, 'trimmer_config_numucc.txt')   # shipped in the code tarball


def follow_grid(sel, a, state, g):
    '''check the submitted jobs: finished ones are checked and moved on, failed or lost ones marked failed'''
    queued = None
    for s in sel:
        st = state.data.get(s['name'], {})
        if st.get('status') != 'submitted':
            continue
        task = st['task']
        kind, marker = G.Grid.result(task)
        if kind == 'done':
            out, info = grid_paths(s)
            log('%s: grid job done (%.0f s, %.0f MB peak memory)' % (s['name'], marker.get('wall', 0), marker.get('maxrss_mb', 0)))
            state.update(s['name'], grid=marker)
            finish(s, out, info, state, a)
        elif kind == 'fail':
            logf = os.path.join(task['marker_dir'], task['name'] + '.log')
            state.update(s['name'], status='failed', grid=marker, error='grid job failed at %s (rc %s), log %s' % (marker.get('stage'), marker.get('rc'), logf))
            log('%s: grid job FAILED at %s, log %s' % (s['name'], marker.get('stage'), logf))
        else:
            if queued is None:
                queued = g.queued()
            if G.Grid.cluster(st.get('jobid')) not in queued:
                state.update(s['name'], status='failed', error='grid job %s left the queue without a result (jobsub_fetchlog for its log)' % st.get('jobid'))
                log('%s: grid job %s lost' % (s['name'], st.get('jobid')))


def submit_grid(sel, a, state, g):
    todo = []
    for s in sel:
        st = state.data.get(s['name'], {})
        if st.get('status') in ('done', 'submitted'):
            continue
        if st.get('no_retry'):
            log('%s: failed its check, not resubmitting (%s)' % (s['name'], st.get('error')))
            continue
        if st.get('tries', 0) >= a.max_tries:
            log('%s: %d submissions failed, not resubmitting (--reset it after a look at the log)' % (s['name'], st['tries']))
            continue
        if os.path.exists(s['pnfs']) or os.path.exists(s['exp']):
            log('%s: an output exists that this script did not finish, not touching it: %s' % (s['name'], s['pnfs']))
            continue
        todo.append(s)
    if not todo:
        return
    tasks = [grid_task(s) for s in todo]
    jobid = g.submit(tasks, 'trim_numucc', **GRID_RES)
    for s, t in zip(todo, tasks):
        if a.dry_run:
            log('%s: %s' % (s['name'], t['cmd'][2]))
            continue
        st = state.data.get(s['name'], {})
        state.update(s['name'], status='submitted', mode='grid', jobid=jobid, task=t, tries=st.get('tries', 0) + 1,
                     error=None, submitted=time.strftime('%Y-%m-%d %H:%M:%S'))


# ---------------------------------------------------------------------------------------------------------------------
def delete_old(sel, state, yes):
    todo = []
    for s in sel:
        if state.data.get(s['name'], {}).get('status') != 'done':
            continue
        if not (os.path.exists(s['pnfs']) and os.path.exists(s['exp']) and os.path.getsize(s['pnfs']) == os.path.getsize(s['exp'])):
            log('%s: the new file is not in place on /pnfs and /exp, keeping the old ones' % s['name'])
            continue
        todo += [p for p in s['old'] if os.path.exists(p)]
    if not todo:
        print('no old checkouts to remove')
        return
    print('old checkouts replaced by their _v2 version (%d files, %.1f GB):' % (len(todo), sum(os.path.getsize(p) for p in todo) / 1e9))
    for p in todo:
        print('  rm %s' % p)
    if not yes and input('remove them? type "yes": ').strip() != 'yes':
        print('nothing removed')
        return
    for p in todo:
        remove(p)
        log('removed %s' % p)


def status(sel, state):
    for s in sel:
        st = state.data.get(s['name'], {})
        c = st.get('check', {})
        extra = ('%d -> %d events' % (c['events_in'], c['events_out'])) if c else ''
        print('%-10s %-6s %-95s %s%s' % (st.get('status', 'pending'), st.get('mode', ''), s['name'], extra,
                                          ('  ' + st['error']) if st.get('error') else ''))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-j', '--jobs', type=int, default=2, help='local: files at the same time (default 2)')
    ap.add_argument('--grid', action='store_true', help='run the trims as grid jobs (and follow the submitted ones)')
    ap.add_argument('--wait', action='store_true', help='--grid: follow the jobs until all are done or failed')
    ap.add_argument('--poll', type=int, default=300, help='--wait: seconds between checks (default 300)')
    ap.add_argument('--max-tries', type=int, default=2, help='--grid: submissions per file before giving up (default 2)')
    ap.add_argument('--match', default='', help='only files whose output name contains this')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--reset', default=None, help='regular expression: forget the progress of the matching files, then exit')
    ap.add_argument('--delete-old', action='store_true', help='remove the old checkouts of the files that are done (asks)')
    ap.add_argument('--yes', action='store_true', help='--delete-old without asking')
    ap.add_argument('--max-loss', type=float, default=0.0, help='fraction of the events, T_pot / SubRun entries and POT the trim may '
                    'remove (default 0: the output must match the input exactly)')
    ap.add_argument('--bin-dir', default=os.path.join(rp.TOP, 'bin'), help='local: the apps')
    ap.add_argument('--app-dir', default=rp.TOP, help='--grid: directory with bin/, lib/, lib64/ for the code tarball')
    ap.add_argument('--config', default=GRID_CONFIG, help='local: the tree_trimmer config')
    ap.add_argument('--state-dir', default=os.path.join(os.getcwd(), 'state_trim_numucc'))
    a = ap.parse_args()
    a.state_dir = os.path.abspath(a.state_dir)
    for d in ('logs', 'info', 'grid'):
        os.makedirs(os.path.join(a.state_dir, d), exist_ok=True)
    state = State(os.path.join(a.state_dir, 'state.json'), dry=a.dry_run or a.status)
    sel = [s for s in samples() if a.match in s['name']]

    if a.status:
        status(sel, state)
        return 0
    lockf = open(os.path.join(a.state_dir, '.lock'), 'w')
    try:
        fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print('another run_trim_numucc.py is using %s' % a.state_dir)
        return 1
    if a.reset:
        names = [n for n in list(state.data) if re.search(a.reset, n)]
        for n in names:
            state.data.pop(n)
        state.save()
        print('forgot %d: %s' % (len(names), ' '.join(names)))
        return 0
    if a.delete_old:
        delete_old(sel, state, a.yes)
        return 0
    missing = [s['input'] for s in sel if not os.path.exists(s['input'])]
    if missing:
        print('missing inputs (run on a machine with /pnfs):\n  ' + '\n  '.join(missing))
        if not a.dry_run:
            return 1

    if a.grid:
        g = G.Grid(os.path.join(a.state_dir, 'grid'), a.app_dir, dry=a.dry_run, log=log)
        while True:
            follow_grid(sel, a, state, g)
            submit_grid(sel, a, state, g)
            pending = [s for s in sel if state.data.get(s['name'], {}).get('status') == 'submitted']
            if not a.wait or a.dry_run or not pending:
                break
            log('%d grid jobs running, next check in %d s' % (len(pending), a.poll))
            time.sleep(a.poll)
        status(sel, state)
        return 0

    # local
    todo = [s for s in sel if state.data.get(s['name'], {}).get('status') != 'done']
    for s in list(todo):
        st = state.data.get(s['name'], {})
        if st.get('status') == 'submitted':
            log('%s: submitted to the grid (%s), follow it with --grid' % (s['name'], st.get('jobid')))
            todo.remove(s)
        elif os.path.exists(s['pnfs']) or os.path.exists(s['exp']):
            log('%s: an output exists that this script did not finish, not touching it' % s['name'])
            todo.remove(s)
    if a.dry_run:
        for s in todo:
            print('%s\n  %s %s %s %s %s\n  -> %s, %s' % (s['name'], os.path.join(a.bin_dir, 'tree_trimmer'), s['input'],
                                                          tmp_name(s['pnfs']), a.config, ' '.join(FLAGS), s['pnfs'], s['exp']))
        print('%d of %d files to do' % (len(todo), len(sel)))
        return 0

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
    with ThreadPoolExecutor(max(1, a.jobs)) as ex:
        ok = list(ex.map(lambda s: run_local(s, a, state), todo))
    log('%d of %d done' % (sum(ok), len(todo)))
    status(sel, state)
    return 0 if all(ok) else 1


if __name__ == '__main__':
    sys.exit(main())
