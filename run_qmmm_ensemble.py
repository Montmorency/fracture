#!/usr/bin/env python
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)

import os
import sys
import subprocess
import socket
import glob
import shutil

from bgqtools import (get_bootable_blocks, boot_blocks, block_corner_iter,
                      set_unbuffered_stdout)

set_unbuffered_stdout()

acct  = 'SiO2_Fracture'
time  = 60
queue = 'default'
mapping = 'ABCDET'
scratch = os.getcwd()
vasp    = '/projects/SiO2_Fracture/iron/vasp.bgq'
envargs = '--envs RUNJOB_MAPPING=%s --envs MPIRUN_ENABLE_TTY_REPORTING=0' % mapping
npj     = 64 # nodes per job
ppn     = 4 # MPI tasks per node

hostname = socket.gethostname()
print 'Hostname: %s' % hostname

jobdirs = glob.glob('T3*')[:2]
print 'jobdirs = %s' % jobdirs

njobs = len(jobdirs)
nodes = npj*njobs

if 'COBALT_PARTSIZE' not in os.environ:
    print 'Not running under control of cobalt. Launching qsub...'
    qsub_args = 'qsub -A %s -n %d -t %d -q %s --mode script --disable_preboot %s' % (acct, nodes, time, queue, ' '.join(sys.argv))
    print qsub_args
    os.system(qsub_args)
    sys.exit(1)

partsize  = int(os.environ['COBALT_PARTSIZE'])
partition = os.environ['COBALT_PARTNAME']
jobid     = int(os.environ['COBALT_JOBID'])

assert nodes == partsize

print 'Nodes per job: %d' % npj
print 'MPI tasks per node: %d' % ppn
print 'Number of sub-block jobs: %d' % njobs
print 'Number of nodes: %d' % nodes

blocks = get_bootable_blocks(partition, nodes)
print 'Available blocks: %s' % blocks

boot_blocks(blocks)

# start sub-block jobs with background runjob helper processes
jobs = []
logs = []
print jobdirs
for job, (block, corner, shape) in zip(jobdirs, block_corner_iter(blocks, npj)):
    print 'JOBDIR and BCS', job, (block, corner, shape)
    shutil.copy(os.path.join(scratch, 'run_qmmm_e.py'), os.path.join(scratch,job))
    os.chdir(os.path.join(scratch, job))
    log = open('%d.vasp.stdout' % jobid, 'w')
# We pass lots of arguments to the run_qmmm script that tell it what to do
# sim temperature, check force error etc.
# and where (i.e. which block corner shape) to do it.
    physargs    = '-st {0}'.format(300)
    locargs     = '--block {0} --corner {1} --shape {2}'.format(block, corner, shape)
    nodeargs    = '--npj {0} --ppn {1}'.format(npj, ppn)
    runjob_args = ('python run_qmmm_e.py {physargs} {locargs} {nodeargs} '.format(physargs=physargs, locargs=locargs, nodeargs=nodeargs)).split()
    jobs.append(subprocess.Popen(runjob_args, stdout=log))
    logs.append(log)
# wait for all background jobs to finish, then flush their logs
for (job, log) in zip(jobs, logs):
    job.wait()
    log.flush()
