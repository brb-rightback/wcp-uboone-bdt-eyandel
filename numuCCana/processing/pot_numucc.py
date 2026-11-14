#!/usr/bin/env python3
# POT and triggers of the numuCC data checkouts, checked against the per-file POT/triggers column of the sample wiki.
#
# Files: the 13 BNB beam-on checkout_<name>_v2.root of run_trim_numucc.py and the BNB beam-off (EXT) checkouts of
# run_production.py, in /exp/uboone/data/users/bbogart/numucc/. For each, numuCCana/pot_info/trimmed_check/trimmed_gdi.py
# makes the run-subrun lists of nuselection/SubRun and wcpselection/T_pot and runs getDataInfo.py -v3 on them (beam on:
# POT tor860_wcut, triggers E1DCNT_wcut; beam off: triggers EXT, with --prescale --nthreads 1).
#
# Wiki check (SubRun list, the one the wiki column uses): each file is compared with the wiki entry it corresponds to:
#   beam on   checkout_<name>_v2.root   <->  the _v2 file it was trimmed from (trim_checkout_<name>_v2.root): the trim
#             kept every subrun, so the POT must agree to the 4 digits the wiki prints and the triggers exactly
#   beam off  checkout_<name>.root      <->  the uboonepro checkout of the same name (wc processed; runs 1, 3: the file
#             without the WC BDT training subruns, as run_production.py removes them too): reported, a difference
#             means the two productions kept different subruns
# Files without a wiki entry (or an entry without POT/triggers) are listed. Exit code 1 if a beam-on file disagrees.
#
# Run where getDataInfo.py works (python2, or the SL7 container wrapper trimmed_gdi.py uses) and uproot or PyROOT reads
# the files:
#   python3 numuCCana/processing/pot_numucc.py [-j 4] [--outdir DIR] [--wiki numuCCana/mcc910wiki.txt]
# Output in --outdir (default ./pot_numucc): trimmed_gdi/ (its lists, getDataInfo.py cache, summary.txt, results.json),
# wiki_check.txt (the comparison), pot_numucc.json ({file: {pot, triggers, subruns, T_pot list values, wiki values}}) for
# the config generators.
import argparse, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_production as rp
import run_trim_numucc as rt

TRIMMED_GDI = os.path.join(rp.TOP, 'numuCCana', 'pot_info', 'trimmed_check', 'trimmed_gdi.py')
POT_REL = 5e-4   # the wiki prints the POT with 4 significant digits


def files():
    '''[(file, kind, wiki file name it is compared with, its /pnfs path for beam on)]'''
    out = [(s['exp'], 'beam on', os.path.basename(s['input']), s['input']) for s in rt.samples()]
    out += [(os.path.join(rt.EXP_OUT, 'checkout_%s.root' % s['name']), 'beam off', 'checkout_%s.root' % s['name'], None)
            for s in rp.SAMPLES if s['kind'] == 'beam_off']
    return out


def read_wiki_values(path):
    '''{file name: (directory, value cell)} of every file entry of the wiki with a POT/triggers cell: a file line is
    "<name>.root | type | value |" under its directory line, or, for an entry with an /exp alternative ("or"), the /pnfs
    file line has no cells and takes the cells of the /exp file line. Only the /pnfs directories are kept.'''
    out, cur, pending = {}, None, []
    for line in open(path).read().split('\n'):
        l = line.strip().strip('|').strip()
        m = re.match(r'(/(?:pnfs|exp)/\S*)', l)
        if m and not l.split('|')[0].strip().endswith('.root'):
            cur = re.sub('/+', '/', m.group(1).rstrip('/'))
            continue
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        tok = cells[0]
        if pending and line.strip().startswith('|') and not tok.endswith('.root') and not re.match(r'^([\\/]\d+\.|_\.)', tok):
            cells, tok = [''] + cells, ''   # the cells of the file above, on a line of their own ("| type | value |")
        elif not tok.endswith('.root') or cur is None:
            continue
        if tok and cur.startswith('/pnfs'):
            pending.append((tok, cur))
        if len(cells) >= 3:
            for name, d in pending:
                out.setdefault(name, []).append((d, cells[2]))
            pending = []
    return out


def parse_value(cell):
    pot = re.search(r'([0-9.]+e[+-]?[0-9]+) POT \((\w+)\)', cell)
    trig = re.search(r'([0-9.]+) triggers \((\w+)\)', cell)
    return (float(pot.group(1)) if pot else None, float(trig.group(1)) if trig else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default=os.path.join(os.getcwd(), 'pot_numucc'))
    ap.add_argument('--wiki', default=os.path.join(rp.TOP, 'numuCCana', 'mcc910wiki.txt'))
    ap.add_argument('-j', '--jobs', type=int, default=4)
    ap.add_argument('--gdi-cmd', default=None, help='passed to trimmed_gdi.py (default: its own)')
    ap.add_argument('--skip-gdi', action='store_true', help='only redo the wiki check from an existing results.json')
    a = ap.parse_args()
    F = files()
    missing = [f for f, _, _, _ in F if not os.path.exists(f)]
    if missing:
        sys.exit('missing files:\n  ' + '\n  '.join(missing))
    gdir = os.path.join(a.outdir, 'trimmed_gdi')
    if not a.skip_gdi:
        cmd = [sys.executable, TRIMMED_GDI, '--outdir', gdir, '-j', str(a.jobs)] + (['--gdi-cmd', a.gdi_cmd] if a.gdi_cmd else []) + [f for f, _, _, _ in F]
        print(' '.join(cmd))
        if subprocess.call(cmd) != 0:
            sys.exit('trimmed_gdi.py failed')
    res = dict((r['path'], r) for r in json.load(open(os.path.join(gdir, 'results.json'))))
    wiki = read_wiki_values(a.wiki)

    lines = ['POT / triggers (getDataInfo.py -v3, SubRun list) vs the per-file column of %s' % a.wiki, '',
             '%-8s %-92s %12s %12s %12s %12s  %s' % ('kind', 'file', 'POT', 'wiki POT', 'triggers', 'wiki trig', 'result')]
    table, bad_on, notes = {}, 0, []
    for f, kind, wname, wpath in F:
        r = res.get(f)
        if r is None:
            lines.append('%-8s %-92s  no getDataInfo.py result' % (kind, os.path.basename(f)))
            bad_on += kind == 'beam on'
            continue
        x, t = r['lists']['SubRun'], r['lists'].get('T_pot', {})
        entry = wiki.get(wname, [])
        if kind == 'beam on':    # the _v2 file itself (its Wire-Cell only twin has the same name)
            entry = [e for e in entry if e[0] == os.path.dirname(re.sub('/+', '/', wpath))]
        else:                    # the wc processed file (not the Wire-Cell only one, which has a _v2 name)
            entry = [e for e in entry if '/wc_only' not in e[0]]
        wpot, wtrig = parse_value(entry[0][1]) if len(entry) == 1 else (None, None)
        if len(entry) != 1:
            result = 'no wiki entry' if not entry else 'several wiki entries (%d)' % len(entry)
        else:
            ok_pot = (x['pot'] is None and wpot is None) or (wpot is not None and x['pot'] is not None and abs(x['pot'] - wpot) <= POT_REL * wpot)
            ok_trig = wtrig is not None and x['trig'] is not None and abs(x['trig'] - wtrig) < 0.5
            result = 'match' if ok_pot and ok_trig else 'DIFFERENT (%s)' % ', '.join(n for n, ok in (('POT', ok_pot), ('triggers', ok_trig)) if not ok)
        if result != 'match' and kind == 'beam on':
            bad_on += 1
        if result != 'match' and kind == 'beam off':
            notes.append(os.path.basename(f))
        lines.append('%-8s %-92s %12s %12s %12s %12s  %s' % (
            kind, os.path.basename(f)[:92], '%.4e' % x['pot'] if x['pot'] else '-', '%.4e' % wpot if wpot else '-',
            '%d' % round(x['trig']) if x['trig'] is not None else '-', '%d' % wtrig if wtrig else '-', result))
        table[f] = dict(kind=kind, pot=x['pot'], triggers=x['trig'], subruns=x['nsub'], events=r['nev'],
                        tpot_list=dict(pot=t.get('pot'), triggers=t.get('trig'), subruns=t.get('nsub')),
                        wiki_file=wname, wiki_pot=wpot, wiki_triggers=wtrig, wiki_check=result)
    lines += ['', 'beam on: %d of %d disagree with the wiki' % (bad_on, sum(1 for _, k, _, _ in F if k == 'beam on')),
              'beam off: %d of %d differ from the wiki entry of the same name%s' % (
                  len(notes), sum(1 for _, k, _, _ in F if k == 'beam off'),
                  ' (different subruns kept by the two productions: look at the lists in trimmed_gdi/lists)' if notes else '')]
    txt = '\n'.join(lines) + '\n'
    open(os.path.join(a.outdir, 'wiki_check.txt'), 'w').write(txt)
    json.dump(table, open(os.path.join(a.outdir, 'pot_numucc.json'), 'w'), indent=1, sort_keys=True)
    print(txt)
    print('wrote %s, %s' % (os.path.join(a.outdir, 'wiki_check.txt'), os.path.join(a.outdir, 'pot_numucc.json')))
    return 1 if bad_on else 0


if __name__ == '__main__':
    sys.exit(main())
