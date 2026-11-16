#!/usr/bin/env python3
# Worker side of grid.py: runs one task of a submitted task list (python 3.6, the SL7 container's python).
#
#   python3 grid_task.py <task list json> <index>          (called by grid_job.sh with $PROCESS)
#
# A task (see grid.py) is a command with placeholders:
#   {in0}, {in1}, ...   its inputs: an xrootd URL (mode "stream", read by the app with TFile::Open) or a local copy in
#                       ./in/ (mode "copy", made with ifdh cp before the command)
#   {out0}, {out1}, ... its outputs: files ./out/<name> the command writes, copied with ifdh cp to their "dest" on /pnfs
#                       (scratch) after the command succeeded
#   {code}              the unpacked code tarball (bin/, lib/, lib64/, processing/)
#   {work}              the working directory of the task (where its "files", {name: content}, are written first)
# The command's stdout/stderr go to <task name>.log, which is copied next to the first output's dest, together with the
# marker <task name>.done (all outputs copied) or <task name>.fail (anything failed): a json with the exit code, the
# wall time, the peak memory of the command and the disk used, which grid.py reads to follow the jobs and to measure
# the resources of a pilot submission. The marker is copied last, so a .done marker means the outputs are complete.
import json, os, resource, shutil, socket, subprocess, sys, time

XROOTD = 'root://fndcadoor.fnal.gov:1094'


def xrootd_url(path):
    '''/pnfs/uboone/... -> root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/...'''
    p = os.path.normpath(path)
    if not p.startswith('/pnfs/'):
        raise ValueError('only /pnfs files can be read on the grid: %s' % path)
    return XROOTD + '/pnfs/fnal.gov/usr/' + p[len('/pnfs/'):]


def sh(cmd, log):
    log.write('+ %s\n' % ' '.join(cmd))
    log.flush()
    return subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)


def du(path):
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


def main():
    tasks = json.load(open(sys.argv[1]))
    task = tasks[int(sys.argv[2])]
    name, code = task['name'], os.environ['CODE_DIR']
    work = os.path.abspath('task_work')
    for d in ('in', 'out'):
        os.makedirs(os.path.join(work, d), exist_ok=True)
    os.chdir(work)
    # bdt_convert reads its BDT weights from ./weights (tarballs made with weights only)
    if os.path.isdir(os.path.join(code, 'processing', 'weights')) and not os.path.exists('weights'):
        os.symlink(os.path.join(code, 'processing', 'weights'), 'weights')
    logname = name + '.log'
    log = open(logname, 'w')
    # cluster: which submission wrote the marker (a job of an earlier submission may still run and write one later)
    info = dict(name=name, host=socket.gethostname(), cluster=os.environ.get('CLUSTER'), start=time.time(), rc=None, stage=None)

    def finish(stage, rc):
        info.update(stage=stage, rc=rc, end=time.time(), wall=time.time() - info['start'],
                    maxrss_mb=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024., disk_mb=du(work) / 1e6)
        ok = stage == 'done'
        marker = name + ('.done' if ok else '.fail')
        json.dump(info, open(marker, 'w'), indent=1)
        log.write('\n%s: %s (rc %s)\n' % (name, stage, rc))
        log.close()
        dest_dir = task['marker_dir']
        for f in (logname, marker):
            subprocess.call(['ifdh', 'rm', os.path.join(dest_dir, f)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.call(['ifdh', 'cp', f, os.path.join(dest_dir, f)])
        sys.exit(0 if ok else 1)

    # small files the task carries (e.g. the results of earlier tasks its check needs), written before the command
    for fname, content in task.get('files', {}).items():
        if os.path.dirname(fname):
            os.makedirs(os.path.join(work, os.path.dirname(fname)), exist_ok=True)
        with open(os.path.join(work, fname), 'w') as f:
            f.write(content)
    subst = {'code': code, 'work': work}
    for i, inp in enumerate(task.get('inputs', [])):
        if inp.get('mode', 'stream') == 'copy':
            local = os.path.join(work, 'in', os.path.basename(inp['path']))
            if sh(['ifdh', 'cp', inp['path'], local], log) != 0 or not os.path.exists(local):
                finish('copy-in %s' % inp['path'], 1)
            subst['in%d' % i] = local
        else:
            subst['in%d' % i] = xrootd_url(inp['path'])
    for i, out in enumerate(task.get('outputs', [])):
        subst['out%d' % i] = os.path.join(work, 'out', out['name'])
    try:
        cmd = [c.format(**subst) for c in task['cmd']]
    except (KeyError, IndexError, ValueError) as e:   # a { } in the command that is not a placeholder
        log.write('bad placeholder in the command: %r\n' % e)
        finish('placeholders', 1)

    t0 = time.time()
    rc = sh(cmd, log)
    info['cmd_wall'] = time.time() - t0
    if rc != 0:
        finish('command', rc)
    for i, out in enumerate(task.get('outputs', [])):
        local = subst['out%d' % i]
        if not os.path.exists(local):
            finish('missing output %s' % out['name'], 1)
        subprocess.call(['ifdh', 'rm', out['dest']], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if sh(['ifdh', 'cp', local, out['dest']], log) != 0:
            finish('copy-out %s' % out['name'], 1)
        info.setdefault('outputs', {})[out['dest']] = os.path.getsize(local)
    finish('done', 0)


if __name__ == '__main__':
    main()
