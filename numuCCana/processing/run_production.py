#!/usr/bin/env python3
# Production of the numuCC BNB checkout (CV) and flux / cross section systematics (XsFlux) files from the unprocessed
# SURPRISE ntuples listed in numuCCana/mcc910wiki.txt. Per sample:
#
#   bdt_convert (per input file) -> [hadd, if the sample has several input files]
#     nu overlay only: prune_weightsep24_trees (all knobs) -> merge_xf ... all -> tree_trimmer (trimmer_config_numucc_xf.txt)
#   convert_cv -> tree_trimmer (trimmer_config_numucc.txt) -> checks -> delete the intermediate files
#
# The run 1 and 3 open data beam-on samples have no unprocessed file of their own: they are made from the new run 1 / 3
# beam-on checkout with tree_trimmer -i full_opendata{1,3}_anti_list.txt (as the production team made them).
# bdt_convert flags follow the production table (numuCCana/proscessingwiki.txt, Production Details), with config.txt
# everywhere (no spline weights): -p1 -r1 -a<samdef> -t<config>, plus -l<trainlist> -g<file type> -b0 for the run 1 / 3 nu
# overlay and dirt (WC BDT training subruns) and -i<zero lifetime list> for the run 4b / 5 nu overlay and dirt.
#
# The configs and lists it uses (config.txt, trimmer_config_numucc[_xf].txt, trainlist_all.dat, zero_e_lieftime_list.txt,
# full_opendata{1,3}_anti_list.txt, weights/) are read from this directory; the apps from <repo>/bin, the samples from
# numuCCana/mcc910wiki.txt.
#
# Run inside the SL7 container with uboonecode set up (the binaries, hadd and PyROOT are needed):
#   python3 numuCCana/processing/run_production.py --dry-run            # commands of every step, nothing is run
#   python3 numuCCana/processing/run_production.py -j 4                 # run (resumes from the state file)
#   python3 numuCCana/processing/run_production.py --status             # progress of each sample
#   python3 numuCCana/processing/run_production.py --kinds nu_overlay --runs 1,4b -j 2
#   python3 numuCCana/processing/run_production.py --reset <sample>     # forget the progress of a sample (or "all")
#
# Resuming: the progress of every step is in <state dir>/state.json (written atomically, on the local disk rather than
# /pnfs). Each step writes to a temporary file that is renamed once the step has finished and passed its check, so a
# killed job never leaves a file that looks complete. On a restart, finished steps are skipped, and a step is re-run if
# its command changed or if an output a remaining step needs was deleted. Ctrl-C / kill stop the running apps; if the
# script was killed hard (kill -9) and an app is still running, the restart refuses to start until it has finished.
# Logs of every step: <state dir>/logs/<sample>/<step>.log (kept after the intermediate files are deleted).
import argparse, fcntl, json, os, re, signal, socket, subprocess, sys, threading, time, traceback
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.abspath(os.path.join(HERE, '../..'))
OUT_DIR = '/pnfs/uboone/persistent/users/bbogart/numucc'

# ---------------------------------------------------------------------------------------------------------------------
# samples
# kind: nu_overlay, dirt, beam_on, beam_off, opendata (beam on); train: file type of the BDT training subruns to remove
# (bdt_convert -g, with -l<trainlist> -b0); zero: zero electron lifetime runs to remove (bdt_convert -i); files: input
# files when the wiki is wrong; parent / anti: open data made from the checkout of parent with tree_trimmer -i<anti>
R123 = 'MCC9.10_Run123_v10_04_07_20_BNB_'
SAMPLES = [
    dict(name=R123 + 'nu_overlay_surprise_reco2_hist_1', kind='nu_overlay', run='1', train='prodgenie_bnb_nu_overlay_run1'),
    dict(name=R123 + 'nu_overlay_surprise_reco2_hist_2', kind='nu_overlay', run='2'),
    dict(name=R123 + 'nu_overlay_surprise_reco2_hist_3', kind='nu_overlay', run='3', train='prodgenie_bnb_nu_overlay_run3'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4a', kind='nu_overlay', run='4a'),
    dict(name='MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist', kind='nu_overlay', run='4b', zero=True),
    dict(name='MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4c', kind='nu_overlay', run='4c'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4d', kind='nu_overlay', run='4d'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_5', kind='nu_overlay', run='5', zero=True),

    dict(name='MCC9.10_Run123_v10_04_07_23_BNB_dirt_overlay_surprise_reco2_hist_1', kind='dirt', run='1', train='prodgenie_bnb_dirt_overlay_run1'),
    dict(name='MCC9.10_Run123_v10_04_07_23_BNB_dirt_overlay_surprise_reco2_hist_2', kind='dirt', run='2'),
    dict(name='MCC9.10_Run123_v10_04_07_23_BNB_dirt_overlay_surprise_reco2_hist_3', kind='dirt', run='3', train='prodgenie_bnb_dirt_overlay_run3'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_20_BNB_dirt_overlay_retuple_retuple_hist_4a', kind='dirt', run='4a'),
    dict(name='MCC9.10_Run4b_v10_04_07_20_BNB_dirt_overlay_retuple_retuple_hist', kind='dirt', run='4b', zero=True),
    dict(name='MCC9.10_Run4acd5_v10_04_07_20_BNB_dirt_overlay_retuple_retuple_hist_4c', kind='dirt', run='4c'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_20_BNB_dirt_overlay_retuple_retuple_hist_4d', kind='dirt', run='4d'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_20_BNB_dirt_overlay_retuple_retuple_hist_5', kind='dirt', run='5', zero=True),

    dict(name=R123 + 'beam_on_data_surprise_reco2_hist_1', kind='beam_on', run='1'),
    dict(name=R123 + 'beam_on_data_surprise_reco2_hist_2', kind='beam_on', run='2'),
    dict(name=R123 + 'beam_on_data_surprise_reco2_hist_3', kind='beam_on', run='3'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_4a', kind='beam_on', run='4a'),
    dict(name='MCC9.10_Run4b_v10_04_07_20_BNB_beam_on_metapatch_retuple_retuple_hist', kind='beam_on', run='4b'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_4c', kind='beam_on', run='4c'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_4d', kind='beam_on', run='4d'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_5', kind='beam_on', run='5'),

    dict(name=R123 + 'beam_on_data_surprise_reco2_hist_1_5e19opendata', kind='opendata', run='1',
         parent=R123 + 'beam_on_data_surprise_reco2_hist_1', anti='full_opendata1_anti_list.txt'),
    dict(name=R123 + 'beam_on_data_surprise_reco2_hist_3_1e19opendata', kind='opendata', run='3',
         parent=R123 + 'beam_on_data_surprise_reco2_hist_3', anti='full_opendata3_anti_list.txt'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_opendata_19550', kind='opendata', run='4a'),
    dict(name='MCC9.10_Run4b_v10_04_07_20_BNB_beam_on_metapatch_retuple_retuple_hist_opendata_20700', kind='opendata', run='4b'),

    dict(name=R123 + 'beam_off_data_surprise_reco2_hist_1', kind='beam_off', run='1'),
    dict(name=R123 + 'beam_off_data_surprise_reco2_hist_2', kind='beam_off', run='2'),
    dict(name=R123 + 'beam_off_data_surprise_reco2_hist_3', kind='beam_off', run='3'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_off_retuple_retuple_hist_4a', kind='beam_off', run='4a'),
    dict(name='MCC9.10_Run4b_v10_04_07_20_BNB_beam_off_metapatch_retuple_retuple_hist', kind='beam_off', run='4b'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_off_retuple_retuple_hist_4c', kind='beam_off', run='4c'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_off_retuple_retuple_hist_4d', kind='beam_off', run='4d'),
    dict(name='MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_off_retuple_retuple_hist_5', kind='beam_off', run='5'),
]
MC_KINDS = ('nu_overlay', 'dirt')

# branches every final file must have (tree: branches); the event trees must all have the same number of entries
EVENT_TREES = {
    'wcpselection/T_eval': ['run', 'subrun', 'event', 'samdef', 'match_found', 'stm_lowenergy'],
    'wcpselection/T_BDTvars': ['numu_cc_flag', 'numu_score', 'nue_score', 'all_veto_score', 'VtxAct_bdt_score'],
    'wcpselection/T_PFeval': ['run', 'reco_nuvtxX', 'mcs_emu_MCS'],
    'wcpselection/T_KINEvars': ['kine_reco_Enu', 'kine_energy_particle'],
    'wcpselection/T_spacepoints': ['Trecchargeblob_spacepoints_q'],
    'nuselection/NeutrinoSelectionFilter': ['run', 'sub', 'evt', 'n_pfps'],
    'lantern/EventTree': ['run', 'subrun', 'event', 'haveReco'],
}
MC_EVAL_BRANCHES = ['weight_cv', 'weight_spline', 'truth_nuEnergy']
POT_TREE = ('wcpselection/T_pot', ['runNo', 'subRunNo', 'pot_tor875'])
XF_KNOBS = ['expskin_FluxUnisim', 'horncurrent_FluxUnisim', 'kminus_PrimaryHadronNormalization',
            'kplus_PrimaryHadronFeynmanScaling', 'kzero_PrimaryHadronSanfordWang', 'nucleoninexsec_FluxUnisim',
            'nucleonqexsec_FluxUnisim', 'nucleontotxsec_FluxUnisim', 'piminus_PrimaryHadronSWCentralSplineVariation',
            'pioninexsec_FluxUnisim', 'pionqexsec_FluxUnisim', 'piontotxsec_FluxUnisim',
            'piplus_PrimaryHadronSWCentralSplineVariation', 'reinteractions_piminus_Geant4',
            'reinteractions_piplus_Geant4', 'reinteractions_proton_Geant4', 'All_UBGenie']
XF_TREE = ('wcpselection/T_weight', ['weight_cv', 'mcweight_filled'] + XF_KNOBS)


class StepError(Exception):
    pass


# ---------------------------------------------------------------------------------------------------------------------
# wiki
def read_wiki(path):
    '''{tuple of input files: tag} of the BNB "unprocessed" table cells of the MCC9.10 sample wiki'''
    txt = open(path).read()
    if 'h2. NuMI Overlay' in txt:
        txt = txt[:txt.index('h2. NuMI Overlay')]
    cells = txt.split('|')
    found = {}
    for i, c in enumerate(cells):
        pairs = re.findall(r'(/pnfs/\S+?/)\s*\n\s*(\S+?\.root)', c)
        if not pairs:
            continue
        j = i + 1
        while j < len(cells) and not cells[j].strip():
            j += 1
        tag = cells[j].strip() if j < len(cells) else ''
        if 'unprocessed' in tag and 'should use' not in tag:
            found[tuple(d + f for d, f in pairs)] = tag
    return found


def wiki_files(sample, wiki):
    '''input files of a sample: a cell with exactly <name>.root, or several files all named <name>_<piece>.root'''
    name = sample['name']
    matches = []
    for files in wiki:
        bases = [os.path.basename(f)[:-len('.root')] for f in files]
        if (len(bases) == 1 and bases[0] == name) or (len(bases) > 1 and all(b.startswith(name + '_') for b in bases)):
            matches.append(list(files))
    if len(matches) != 1:
        raise StepError('%d unprocessed entries for %s in the wiki: %s' % (len(matches), name, matches))
    return matches[0]


# ---------------------------------------------------------------------------------------------------------------------
# ROOT checks (PyROOT is imported only when needed, and used by one thread at a time)
ROOT_LOCK = threading.Lock()
_ROOT = None


def root():
    global _ROOT
    if _ROOT is None:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        ROOT.gErrorIgnoreLevel = ROOT.kError
        _ROOT = ROOT
    return _ROOT


def inspect_file(path, trees, pot=False):
    '''{tree: entries} of the trees (each must exist with the branches listed), plus the summed pot_tor875 of T_pot'''
    with ROOT_LOCK:
        R = root()
        f = R.TFile.Open(path)
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
            if pot:
                t = f.Get(POT_TREE[0])
                n = int(t.GetEntries())
                total = 0.
                if n > 0:
                    t.SetEstimate(n + 1)
                    t.Draw('pot_tor875', '', 'goff')
                    v = t.GetV1()
                    total = sum(v[i] for i in range(n))
                info['pot'] = total
        finally:
            f.Close()
        return info


def check_events(info, path):
    n = [info[t] for t in EVENT_TREES]
    if n[0] <= 0:
        raise StepError('%s: no events' % path)
    if len(set(n)) != 1:
        raise StepError('%s: event trees with different entries %s' % (path, dict((t, info[t]) for t in EVENT_TREES)))


def final_trees(sample, xf=False):
    trees = dict(EVENT_TREES)
    if sample['kind'] in MC_KINDS:
        trees['wcpselection/T_eval'] = EVENT_TREES['wcpselection/T_eval'] + MC_EVAL_BRANCHES
    trees[POT_TREE[0]] = POT_TREE[1]
    if xf:
        trees[XF_TREE[0]] = XF_TREE[1]
    return trees


# ---------------------------------------------------------------------------------------------------------------------
# steps of a sample: name, command (None for steps done in python), inputs, outputs, check
def build_steps(sample, wiki, a):
    name, kind = sample['name'], sample['kind']
    mc = kind in MC_KINDS
    work = os.path.join(a.work_dir, name)
    app = lambda x: os.path.join(a.bin_dir, x)
    final_cv = os.path.join(a.out_dir, 'checkout_%s.root' % name)
    final_xf = os.path.join(a.xf_dir, 'xf_all_%s.root' % name)
    steps = []

    def add(step, cmd, inputs, outputs, check=None, tmp=True, content=None):
        steps.append(dict(step=step, cmd=cmd, inputs=inputs, outputs=outputs, check=check, tmp=tmp, content=content))

    def cv_check(path, out, upstream=None, final=False):
        info = inspect_file(path, final_trees(sample) if final else dict(list(EVENT_TREES.items())[:1] + [POT_TREE]), pot=True)
        if final:
            check_events(info, path)
        elif info['wcpselection/T_eval'] <= 0:
            raise StepError('%s: no events' % path)
        if info[POT_TREE[0]] <= 0:
            raise StepError('%s: empty T_pot' % path)
        if mc and info['pot'] <= 0:
            raise StepError('%s: POT %g' % (path, info['pot']))
        if upstream is not None:
            up = out[upstream]['info']
            if info['pot'] > up['pot'] * (1 + 1e-6):
                raise StepError('%s: POT %g above the input POT %g' % (path, info['pot'], up['pot']))
            if info['wcpselection/T_eval'] > up['wcpselection/T_eval']:
                raise StepError('%s: more events than the input' % path)
        return info

    if sample.get('parent'):
        parent_cv = os.path.join(a.out_dir, 'checkout_%s.root' % sample['parent'])
        add('trim_cv', [app('tree_trimmer'), parent_cv, '{tmp}', a.trim_config, '-o1', '-k1', '-c1',
                        '-i' + os.path.join(HERE, sample['anti'])], [parent_cv], [final_cv],
            lambda p, out: cv_check(p, out, final=True))
        add('cleanup', None, [], [], tmp=False)
        return steps

    files = sample.get('files') or wiki_files(sample, wiki)
    flags = ['-p1', '-r1', '-a' + name, '-t' + a.config]
    if sample.get('train'):
        flags += ['-l' + a.trainlist, '-g' + sample['train'], '-b0']
    if sample.get('zero'):
        flags += ['-i' + a.zero_list]
    bdt = os.path.join(work, 'bdt_convert_%s.root' % name)
    pieces, piece_steps = [], []
    for f in files:
        piece = os.path.basename(f)[len(name):-len('.root')]  # '' for a single input file, e.g. '_02' otherwise
        out = bdt if len(files) == 1 else os.path.join(work, 'bdt_convert_%s%s.root' % (name, piece))
        add('bdt_convert' + piece, [app('bdt_convert'), f, '{tmp}'] + flags, [f], [out], lambda p, out: cv_check(p, out))
        pieces.append(out)
        piece_steps.append('bdt_convert' + piece)
    if len(files) > 1:
        def hadd_check(p, out):
            info = cv_check(p, out)
            n = sum(out[x]['info']['wcpselection/T_eval'] for x in piece_steps)
            if info['wcpselection/T_eval'] != n:
                raise StepError('%s: %d events, %d in the inputs' % (p, info['wcpselection/T_eval'], n))
            return info
        add('hadd', ['hadd', '-f', '{tmp}'] + pieces, pieces, [bdt], hadd_check)

    if kind == 'nu_overlay':
        weights = os.path.join(work, 'weights_%s.root' % name)
        flist = os.path.join(work, 'prune_input_%s.txt' % name)
        add('prune_list', None, [], [flist], tmp=False, content=''.join(f + '\n' for f in files))
        add('prune', [app('prune_weightsep24_trees'), flist, '{tmp}', 'all'], [flist] + files, [weights],
            lambda p, out: inspect_file(p, {'all': ['run', 'subrun', 'event']}))
        xf = os.path.join(work, 'xf_all_%s.root' % name)
        add('merge_xf', [app('merge_xf'), bdt, weights, '{tmp}', 'all', '-t' + a.config], [bdt, weights], [xf],
            lambda p, out: dict(cv_check(p, out), **inspect_file(p, dict([XF_TREE]))))
    checkout = os.path.join(work, 'checkout_%s.root' % name)
    up = 'hadd' if len(files) > 1 else 'bdt_convert'
    add('convert_cv', [app('convert_cv'), bdt, '{tmp}', '-t' + a.config], [bdt], [checkout],
        lambda p, out: cv_check(p, out, upstream=up))
    data_flag = '-c0' if mc else '-c1'
    add('trim_cv', [app('tree_trimmer'), checkout, '{tmp}', a.trim_config, '-o1', '-k1', data_flag], [checkout], [final_cv],
        lambda p, out: cv_check(p, out, upstream='convert_cv', final=True))
    if kind == 'nu_overlay':
        def xf_check(p, out):
            info = cv_check(p, out, upstream='merge_xf', final=True)
            info.update(inspect_file(p, dict([XF_TREE])))
            if info[XF_TREE[0]] != info['wcpselection/T_eval']:
                raise StepError('%s: T_weight %d entries, T_eval %d' % (p, info[XF_TREE[0]], info['wcpselection/T_eval']))
            return info
        add('trim_xf', [app('tree_trimmer'), xf, '{tmp}', a.trim_config_xf, '-o1', '-k1', '-c0'], [xf], [final_xf], xf_check)
    add('cleanup', None, [], [], tmp=False)
    return steps


def finals(sample, a):
    out = [os.path.join(a.out_dir, 'checkout_%s.root' % sample['name'])]
    if sample['kind'] == 'nu_overlay':
        out.append(os.path.join(a.xf_dir, 'xf_all_%s.root' % sample['name']))
    return out


def tmp_name(path):
    return os.path.join(os.path.dirname(path), 'tmp_' + os.path.basename(path))


def cmd_string(step):
    if step['cmd'] is None:
        return step['step']
    return ' '.join(x.replace('{tmp}', tmp_name(step['outputs'][0])) for x in step['cmd'])


# ---------------------------------------------------------------------------------------------------------------------
# state
class State:
    def __init__(self, path, dry=False):
        self.path, self.dry = path, dry
        self.lock = threading.Lock()
        self.data = {}
        if os.path.exists(path):
            with open(path) as f:
                self.data = json.load(f)

    def save(self):
        if self.dry:
            return
        tmp = self.path + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(self.data, f, indent=1, sort_keys=True)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, self.path)

    def sample(self, name):
        with self.lock:
            return self.data.setdefault(name, {'status': 'pending', 'steps': {}})

    def update(self, name, step=None, **kw):
        with self.lock:
            s = self.data.setdefault(name, {'status': 'pending', 'steps': {}})
            (s['steps'].setdefault(step, {}) if step else s).update(kw)
            self.save()


def log(msg):
    sys.stdout.write('%s %s\n' % (time.strftime('%Y-%m-%d %H:%M:%S'), msg))
    sys.stdout.flush()


RUNNING = {}  # pid -> Popen, killed on SIGINT / SIGTERM
STOP = threading.Event()


def remove(path):
    if os.path.lexists(path):
        os.remove(path)


def pid_alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def run_sample(sample, wiki, a, state):
    name = sample['name']
    if STOP.is_set():
        return False
    s = state.sample(name)
    if s['status'] == 'done':
        if all(os.path.exists(f) for f in finals(sample, a)):
            log('%s: done' % name)
            return True
        log('%s: final files missing, starting the sample again' % name)
        state.update(name, steps={})
        s = state.sample(name)
    try:
        steps = build_steps(sample, wiki, a)
    except StepError as e:
        state.update(name, status='failed', error=str(e))
        log('%s: FAILED %s' % (name, e))
        return False
    st = s['steps']
    done = lambda x: st.get(x['step'], {}).get('status') == 'done' and st[x['step']].get('cmd') == cmd_string(x)

    # first step to (re)run: the first one not done, moved back while a remaining step needs an output that is gone
    first = next((i for i, x in enumerate(steps) if not done(x)), len(steps))
    while True:
        producer = dict((o, i) for i, x in enumerate(steps) for o in x['outputs'])
        missing = [producer[i] for x in steps[first:] for i in x['inputs'] if i in producer and producer[i] < first and not os.path.exists(i)]
        if not missing:
            break
        first = min(missing)
    for x in steps[first:]:
        if x['step'] in st and st[x['step']].get('status') == 'done':
            state.update(name, x['step'], status='pending')

    os.makedirs(os.path.join(a.state_dir, 'logs', name), exist_ok=True)
    state.update(name, status='running', error=None)
    for x in steps[first:]:
        if STOP.is_set():
            return False
        step = x['step']
        logfile = os.path.join(a.state_dir, 'logs', name, step + '.log')
        t0 = time.time()
        state.update(name, step, status='running', cmd=cmd_string(x), started=time.strftime('%Y-%m-%d %H:%M:%S'),
                     host=socket.gethostname(), pid=None, log=logfile)
        log('%s: %s' % (name, step))
        try:
            for o in x['outputs']:
                os.makedirs(os.path.dirname(o), exist_ok=True)
                if x['tmp']:
                    remove(tmp_name(o))
                remove(o)
            if step == 'prune_list':
                with open(x['outputs'][0], 'w') as f:
                    f.write(x['content'])
                info = {}
            elif step == 'cleanup':
                for o in finals(sample, a):
                    if not os.path.exists(o):
                        raise StepError('final file %s missing, not deleting the intermediate files' % o)
                work = os.path.join(a.work_dir, name)
                if not a.keep_work and os.path.isdir(work):
                    for f in os.listdir(work):
                        remove(os.path.join(work, f))
                    os.rmdir(work)
                info = {}
            else:
                tmp = tmp_name(x['outputs'][0])
                cmd = [c.replace('{tmp}', tmp) for c in x['cmd']]
                with open(logfile, 'w') as lf:
                    lf.write(' '.join(cmd) + '\n\n')
                    lf.flush()
                    p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=a.run_dir, start_new_session=True)
                    RUNNING[p.pid] = p
                    state.update(name, step, pid=p.pid)
                    rc = p.wait()
                    RUNNING.pop(p.pid, None)
                if STOP.is_set():
                    raise StepError('interrupted')
                if rc != 0:
                    raise StepError('exit code %d, see %s' % (rc, logfile))
                if not os.path.exists(tmp) or os.path.getsize(tmp) == 0:
                    raise StepError('no output %s, see %s' % (tmp, logfile))
                out = dict((k, v) for k, v in state.sample(name)['steps'].items())
                info = x['check'](tmp, out) if x['check'] else {}
                os.replace(tmp, x['outputs'][0])
            state.update(name, step, status='done', pid=None, info=info, seconds=round(time.time() - t0),
                         finished=time.strftime('%Y-%m-%d %H:%M:%S'))
        except Exception as e:
            msg = str(e) if isinstance(e, StepError) else traceback.format_exc()
            state.update(name, step, status='failed', pid=None, error=msg)
            state.update(name, status='interrupted' if STOP.is_set() else 'failed', error='%s: %s' % (step, msg))
            log('%s: %s %s: %s' % (name, step, 'stopped' if STOP.is_set() else 'FAILED', msg.strip()))
            return False
    state.update(name, status='done', finished=time.strftime('%Y-%m-%d %H:%M:%S'))
    log('%s: done' % name)
    return True


# ---------------------------------------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__ if __doc__ else 'numuCC production', formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-j', '--jobs', type=int, default=1, help='samples processed at the same time (default 1)')
    p.add_argument('--dry-run', action='store_true', help='print the commands, run nothing')
    p.add_argument('--status', action='store_true', help='print the progress of each sample and exit')
    p.add_argument('--reset', default=None, help='comma separated samples (or "all") whose progress is forgotten, then exit')
    p.add_argument('--samples', default=None, help='regular expression on the sample names')
    p.add_argument('--kinds', default=None, help='comma separated: nu_overlay, dirt, beam_on, beam_off, opendata')
    p.add_argument('--runs', default=None, help='comma separated: 1, 2, 3, 4a, 4b, 4c, 4d, 5')
    p.add_argument('--keep-work', action='store_true', help='do not delete the intermediate files')
    p.add_argument('--out-dir', default=OUT_DIR)
    p.add_argument('--xf-dir', default=os.path.join(OUT_DIR, 'XsFlux'))
    p.add_argument('--work-dir', default=os.path.join(OUT_DIR, 'work'))
    p.add_argument('--state-dir', default=os.path.join(HERE, 'state'), help='state file and logs (keep it off /pnfs)')
    p.add_argument('--bin-dir', default=os.path.join(TOP, 'bin'))
    p.add_argument('--wiki', default=os.path.join(TOP, 'numuCCana/mcc910wiki.txt'))
    p.add_argument('--config', default=os.path.join(HERE, 'config.txt'))
    p.add_argument('--trim-config', default=os.path.join(HERE, 'trimmer_config_numucc.txt'))
    p.add_argument('--trim-config-xf', default=os.path.join(HERE, 'trimmer_config_numucc_xf.txt'))
    p.add_argument('--trainlist', default=os.path.join(HERE, 'trainlist_all.dat'))
    p.add_argument('--zero-list', default=os.path.join(HERE, 'zero_e_lieftime_list.txt'))
    p.add_argument('--weights', default=os.path.join(HERE, 'weights'), help='BDT weights directory of bdt_convert')
    a = p.parse_args()
    a.state_dir = os.path.abspath(a.state_dir)
    a.run_dir = os.path.join(a.state_dir, 'run')

    sel = SAMPLES
    if a.samples:
        sel = [s for s in sel if re.search(a.samples, s['name'])]
    if a.kinds:
        sel = [s for s in sel if s['kind'] in a.kinds.split(',')]
    if a.runs:
        sel = [s for s in sel if s['run'] in a.runs.split(',')]

    os.makedirs(a.state_dir, exist_ok=True)
    state = State(os.path.join(a.state_dir, 'state.json'), dry=a.dry_run or a.status)

    if a.status:
        for s in sel:
            st = state.data.get(s['name'], {})
            steps = st.get('steps', {})
            last = [k for k, v in steps.items() if v.get('status') in ('running', 'failed')]
            print('%-12s %-4s %-12s %s %s' % (s['kind'], s['run'], st.get('status', 'pending'), s['name'],
                                              ('[%s] %s' % (','.join(last), st.get('error') or '')) if st.get('status') != 'done' and last else ''))
        return 0

    lockf = open(os.path.join(a.state_dir, '.lock'), 'w')
    if not a.dry_run:
        try:
            fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            print('another run_production.py is using %s' % a.state_dir)
            return 1

    if a.reset:
        names = [s['name'] for s in SAMPLES] if a.reset == 'all' else a.reset.split(',')
        for n in names:
            state.data.pop(n, None)
        state.save()
        print('reset %d samples' % len(names))
        return 0

    wiki = read_wiki(a.wiki)

    if a.dry_run:
        for s in sel:
            st = state.data.get(s['name'], {}).get('steps', {})
            print('\n=== %s (%s, run %s) %s' % (s['name'], s['kind'], s['run'], state.data.get(s['name'], {}).get('status', 'pending')))
            try:
                steps = build_steps(s, wiki, a)
                made = set(o for y in steps for o in y['outputs'])
                for x in steps:
                    print('  [%-7s] %-18s %s' % (st.get(x['step'], {}).get('status', 'pending'), x['step'], cmd_string(x) if x['cmd'] else ''))
                    for i in x['inputs']:
                        if i not in made and not s.get('parent') and not os.path.exists(i):
                            print('      input not found (on this machine): %s' % i)
            except StepError as e:
                print('  ERROR %s' % e)
        return 0

    # a step still running from a run that was killed hard
    host = socket.gethostname()
    for n, s in state.data.items():
        for k, v in s.get('steps', {}).items():
            if v.get('status') == 'running' and v.get('pid'):
                if v.get('host') == host and pid_alive(v['pid']):
                    print('%s %s is still running (pid %d) from an earlier run, wait for it or kill it' % (n, k, v['pid']))
                    return 1
                if v.get('host') != host:
                    print('WARNING: %s %s was running on %s (pid %d), make sure it is no longer running' % (n, k, v.get('host'), v['pid']))

    # the inputs and tools
    problems = []
    for x in ['bdt_convert', 'convert_cv', 'tree_trimmer', 'merge_xf', 'prune_weightsep24_trees']:
        if not os.access(os.path.join(a.bin_dir, x), os.X_OK):
            problems.append('no %s in %s' % (x, a.bin_dir))
    for f in [a.config, a.trim_config, a.trim_config_xf, a.trainlist, a.zero_list, a.weights] + [os.path.join(HERE, s['anti']) for s in sel if s.get('anti')]:
        if not os.path.exists(f):
            problems.append('missing %s' % f)
    if subprocess.call('command -v hadd > /dev/null', shell=True) != 0:
        problems.append('no hadd (set up uboonecode)')
    try:
        root()
    except Exception as e:
        problems.append('no PyROOT: %s' % e)
    for s in sel:
        if s.get('parent'):
            continue
        try:
            for f in s.get('files') or wiki_files(s, wiki):
                if not os.path.exists(f):
                    problems.append('input of %s not found: %s' % (s['name'], f))
        except StepError as e:
            problems.append(str(e))
    if problems:
        print('\n'.join(problems))
        return 1

    for d in [a.out_dir, a.xf_dir, a.work_dir, a.run_dir]:
        os.makedirs(d, exist_ok=True)
    link = os.path.join(a.run_dir, 'weights')  # bdt_convert reads weights/*.xml from its working directory
    if os.path.realpath(link) != os.path.realpath(a.weights):
        remove(link)
        os.symlink(a.weights, link)

    def stop(signum, frame):
        if not STOP.is_set():
            log('stopping (signal %d), killing the running apps' % signum)
        STOP.set()
        for pid in list(RUNNING):
            try:
                os.killpg(pid, signal.SIGTERM)
            except OSError:
                pass
    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)

    # samples made from other samples go last, after their parent
    first = [s for s in sel if not s.get('parent')]
    second = [s for s in sel if s.get('parent')]
    ok = {}
    with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        for s, r in zip(first, pool.map(lambda s: run_sample(s, wiki, a, state), first)):
            ok[s['name']] = r
        todo = []
        for s in (second if not STOP.is_set() else []):
            parent_ok = ok.get(s['parent'], state.data.get(s['parent'], {}).get('status') == 'done')
            if parent_ok and os.path.exists(os.path.join(a.out_dir, 'checkout_%s.root' % s['parent'])):
                todo.append(s)
            elif not STOP.is_set():
                log('%s: skipped, %s is not done' % (s['name'], s['parent']))
                ok[s['name']] = False
        for s, r in zip(todo, pool.map(lambda s: run_sample(s, wiki, a, state), todo)):
            ok[s['name']] = r

    bad = [n for n, r in ok.items() if not r]
    log('%d of %d samples done%s' % (len(ok) - len(bad), len(ok), (', not done: ' + ', '.join(bad)) if bad else ''))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
