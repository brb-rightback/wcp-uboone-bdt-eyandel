#!/bin/bash
# Grid job of grid.py (one task of a task list per job, $PROCESS = task index), run in the fnal-wn-sl7 container:
#   grid_job.sh <task list file name>
# The code tarball (bin/, lib/, lib64/, processing/) is unpacked by jobsub at $INPUT_TAR_DIR_LOCAL, the task list
# (jobsub -f) is in $CONDOR_DIR_INPUT.
echo "grid_job.sh on $(hostname) at $(date), process ${PROCESS}, task list $1"
export CODE_DIR=${INPUT_TAR_DIR_LOCAL}
if [ ! -d "${CODE_DIR}/bin" ]; then
  echo "no code tarball at ${CODE_DIR}"; exit 2
fi

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v10_04_07_04 -q e26:prof || { echo "setup uboonecode failed"; exit 2; }
export PATH=${CODE_DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${CODE_DIR}/lib:${CODE_DIR}/lib64:${LD_LIBRARY_PATH}

cd ${_CONDOR_SCRATCH_DIR:-$PWD}
python3 ${CODE_DIR}/processing/grid_task.py ${CONDOR_DIR_INPUT}/$1 ${PROCESS}
rc=$?
echo "grid_task.py exit code ${rc} at $(date)"
exit ${rc}
