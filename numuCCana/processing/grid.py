#!/usr/bin/env python3
# Grid (jobsub_lite) execution of tasks for the production and validation scripts, run on a gpvm / build machine with
# /pnfs mounted. A task is one command with its inputs and outputs:
#
#   dict(name='trim_run1',                                      # unique within a submission (file names below)
#        cmd=['{code}/bin/tree_trimmer', '{in0}', '{out0}', ...],  # placeholders: see grid_task.py
#        inputs=[dict(path='/pnfs/...', mode='stream')],         # stream: xrootd URL (TFile::Open); copy: ifdh cp first
#        outputs=[dict(name='x.root', dest='/pnfs/uboone/scratch/users/<user>/.../x.root')],
#        marker_dir='/pnfs/uboone/scratch/...')                  # where <name>.log/.done/.fail go (default: first dest)
#
# Grid.submit() runs a list of tasks as one array job (jobsub_submit -N, task index = $PROCESS) in the SL7 container
# (fnal-wn-sl7), on FermiGrid only (--onsite), with the code tarball (bin/, lib/, lib64/, processing/ scripts and
# configs, optionally processing/weights/) through the RCDS dropbox and the task list as an input file. grid_job.sh
# sets up uboonecode from cvmfs and runs grid_task.py, which runs the task, copies its outputs to their dest (ifdh) and
# writes <marker_dir>/<name>.done or .fail (a json: exit code, failed stage, wall time, peak memory, disk used).
# Grid.result() reads the marker; Grid.queued() lists the clusters still in the queue, so a task with no marker whose
# cluster has left the queue is lost (resubmit it). The calling script checks the outputs and moves them on.
#
# Resources: disk always 30 GB (the slot granularity), memory and expected lifetime per submission.
#
# Self test (submits one small job: reads a /pnfs file over xrootd with ROOT, writes a file back to scratch):
#   python3 numuCCana/processing/grid.py --selftest <a /pnfs root file> [--app-dir DIR]
#   python3 numuCCana/processing/grid.py --selftest-status
import argparse, getpass, glob, hashlib, json, os, re, subprocess, sys, tarfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.abspath(os.path.join(HERE, '../..'))
IMAGE = '/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest'
GROUP = 'uboone'
DISK = '30GB'
SCRATCH = '/pnfs/uboone/scratch/users/%s' % getpass.getuser()
# processing/ files shipped to the workers (scripts, tree configs and lists; not the state directories or logs)
SHIP = ('*.py', '*.sh', '*.txt', '*.dat')


def scratch_dir(*parts):
    return os.path.join(SCRATCH, 'grid', *parts)


class Grid:
    def __init__(self, work_dir, app_dir=TOP, weights=False, dry=False, log=print):
        '''work_dir: local directory for the tarballs and the submitted task lists (e.g. <state dir>/grid)'''
        self.work_dir, self.app_dir, self.weights, self.dry, self.log = work_dir, os.path.abspath(app_dir), weights, dry, log
        os.makedirs(os.path.join(work_dir, 'submissions'), exist_ok=True)
        self._tar = None

    # -------------------------------------------------------------------------------------------- code tarball
    def _members(self):
        '''(path on disk, path in the tarball) of everything shipped'''
        out = []
        for d in ('bin', 'lib', 'lib64'):
            src = os.path.join(self.app_dir, d)
            for root, _, files in os.walk(src, followlinks=True):
                for f in files:
                    p = os.path.join(root, f)
                    out.append((p, os.path.join(d, os.path.relpath(p, src))))
        for pat in SHIP:
            for p in glob.glob(os.path.join(HERE, pat)):
                out.append((p, os.path.join('processing', os.path.basename(p))))
        if self.weights:
            for root, _, files in os.walk(os.path.join(HERE, 'weights'), followlinks=True):
                for f in files:
                    p = os.path.join(root, f)
                    out.append((p, os.path.join('processing', os.path.relpath(p, HERE))))
        if not any(m[1].startswith('bin/') for m in out):
            raise RuntimeError('no apps in %s/bin' % self.app_dir)
        return sorted(out, key=lambda m: m[1])

    def tarball(self):
        '''the code tarball, made again only if a shipped file changed (name = hash of names, sizes and times)'''
        if self._tar:
            return self._tar
        members = self._members()
        h = hashlib.md5()
        for p, arc in members:
            st = os.stat(p)
            h.update(('%s %d %d\n' % (arc, st.st_size, int(st.st_mtime))).encode())
        path = os.path.join(self.work_dir, 'code_%s.tar.gz' % h.hexdigest()[:12])
        if not os.path.exists(path) and not self.dry:
            self.log('grid: making the code tarball %s (%d files)' % (path, len(members)))
            tmp = path + '.tmp'
            with tarfile.open(tmp, 'w:gz') as t:
                for p, arc in members:
                    t.add(p, arcname=arc)
            os.replace(tmp, path)
        self._tar = path
        return path

    # -------------------------------------------------------------------------------------------- submission
    def submit(self, tasks, tag, memory='2GB', lifetime='8h', disk=DISK):
        '''submit tasks as one array job; returns the job id (<cluster>.0@<schedd>), None in a dry run'''
        if not tasks:
            return None
        names = [t['name'] for t in tasks]
        if len(set(names)) != len(names):
            raise RuntimeError('task names must be unique in a submission')
        for t in tasks:
            t.setdefault('marker_dir', os.path.dirname(t['outputs'][0]['dest']) if t.get('outputs') else scratch_dir('markers'))
            for o in t.get('outputs', []):
                if not o['dest'].startswith('/pnfs/'):
                    raise RuntimeError('%s: grid outputs go to /pnfs: %s' % (t['name'], o['dest']))
        if not self.dry:   # ifdh cp does not make directories on /pnfs
            for t in tasks:
                os.makedirs(t['marker_dir'], exist_ok=True)
                for o in t.get('outputs', []):
                    os.makedirs(os.path.dirname(o['dest']), exist_ok=True)
        stamp = time.strftime('%Y%m%d_%H%M%S')
        tfile = os.path.join(self.work_dir, 'submissions', '%s_%s.json' % (tag, stamp))
        with open(tfile, 'w') as f:
            json.dump(tasks, f, indent=1)
        cmd = ['jobsub_submit', '-G', GROUP, '-N', str(len(tasks)), '--onsite', '--singularity-image', IMAGE,
               '--memory', memory, '--disk', disk, '--expected-lifetime', lifetime,
               '--tar_file_name', 'dropbox://' + self.tarball(), '-f', 'dropbox://' + tfile,
               'file://' + os.path.join(HERE, 'grid_job.sh'), os.path.basename(tfile)]
        self.log('grid: %s' % ' '.join(cmd))
        if self.dry:
            return None
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        m = re.search(r'Use job id (\S+) to retrieve output', p.stdout)
        if p.returncode != 0 or not m:
            raise RuntimeError('jobsub_submit failed (%d):\n%s' % (p.returncode, p.stdout[-3000:]))
        jobid = m.group(1)
        with open(tfile + '.jobid', 'w') as f:
            f.write(jobid + '\n')
        self.log('grid: submitted %d tasks as %s' % (len(tasks), jobid))
        return jobid

    # -------------------------------------------------------------------------------------------- following jobs
    @staticmethod
    def result(task):
        '''the marker of a task: ('done' | 'fail', marker json) or (None, None) if there is none yet'''
        for kind in ('done', 'fail'):
            p = os.path.join(task['marker_dir'], '%s.%s' % (task['name'], kind))
            if os.path.exists(p):
                try:
                    return kind, json.load(open(p))
                except ValueError:
                    return None, None   # being written
        return None, None

    def queued(self):
        '''set of the job clusters (<cluster>@<schedd>) of this user still in the queue (idle, running, held)'''
        p = subprocess.run(['jobsub_q', '-G', GROUP, '--user', getpass.getuser()], stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, universal_newlines=True)
        if p.returncode != 0:
            raise RuntimeError('jobsub_q failed:\n%s' % p.stdout[-2000:])
        out = set()
        for m in re.finditer(r'^(\d+)\.\d+@(\S+)', p.stdout, re.M):
            out.add('%s@%s' % (m.group(1), m.group(2)))
        return out

    @staticmethod
    def cluster(jobid):
        m = re.match(r'(\d+)\.\d+@(\S+)', jobid or '')
        return '%s@%s' % (m.group(1), m.group(2)) if m else None

    def fetchlog(self, jobid, dest):
        os.makedirs(dest, exist_ok=True)
        return subprocess.call(['jobsub_fetchlog', '-G', GROUP, '--jobid', jobid, '--destdir', dest])


# ------------------------------------------------------------------------------------------------ self test
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', metavar='PNFS_ROOT_FILE', help='submit one test job reading this file over xrootd')
    ap.add_argument('--selftest-status', action='store_true', help='show the result of the test job')
    ap.add_argument('--app-dir', default=TOP, help='directory with bin/, lib/, lib64/ (default: this repository)')
    ap.add_argument('--work-dir', default=os.path.join(os.getcwd(), 'state_grid_selftest'))
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    g = Grid(a.work_dir, a.app_dir, dry=a.dry_run)
    out = scratch_dir('selftest', 'selftest_copy.root')
    task = dict(name='selftest', marker_dir=scratch_dir('selftest'),
                cmd=['bash', '-c', 'echo apps: $(ls {code}/bin | wc -l); which hadd ifdh; '
                     'root -l -b -q -e "TFile *f = TFile::Open(\\"{in0}\\"); f->ls(); return f->IsZombie();" && '
                     'hadd -f {out0} {in0} && ls -l {out0}'],
                inputs=[dict(path=a.selftest or '/pnfs/none', mode='stream')],
                outputs=[dict(name='selftest_copy.root', dest=out)])
    if a.selftest:
        g.submit([task], 'selftest', memory='2GB', lifetime='1h')
        print('when it has run: python3 %s --selftest-status --work-dir %s' % (__file__, a.work_dir))
    elif a.selftest_status:
        kind, info = Grid.result(task)
        print('marker: %s %s' % (kind, json.dumps(info, indent=1) if info else ''))
        if kind:
            print('log: %s' % os.path.join(task['marker_dir'], 'selftest.log'))
        else:
            jid = sorted(glob.glob(os.path.join(a.work_dir, 'submissions', 'selftest_*.jobid')))
            if jid:
                c = Grid.cluster(open(jid[-1]).read().strip())
                print('job %s %s' % (c, 'in the queue' if c in g.queued() else 'not in the queue (lost? check jobsub_fetchlog)'))
    else:
        ap.print_help()


if __name__ == '__main__':
    main()
