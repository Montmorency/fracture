#!/usr/bin/env python

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

import os
import sys
import subprocess
import socket
import glob

from bgqtools import (get_bootable_blocks, boot_blocks, block_corner_iter,
                      set_unbuffered_stdout)

set_unbuffered_stdout()

acct  = 'SiO2_Fracture'
time  = 720
queue = 'default'
mapping = 'ABCDET'
scratch = os.getcwd()
#vasp = '/home/kermode/exe/vasp.O3.cplx'
vasp  = '/projects/SiO2_Fracture/iron/vasp.bgq'
envargs = '--envs RUNJOB_MAPPING=%s --envs MPIRUN_ENABLE_TTY_REPORTING=0' % mapping
npj = 512 # nodes per job
ppn = 4 # MPI tasks per node

hostname = socket.gethostname()
print 'Hostname: %s' % hostname

#jobdirs = glob.glob('vasp-ref-???')

jobdirs = glob.glob('b*shift*')
print 'jobdirs = %s' % jobdirs

njobs = len(jobdirs)
nodes = npj*njobs

if 'COBALT_PARTSIZE' not in os.environ:
    print 'Not running under control of cobalt. Launching qsub...'
    qsub_args = 'qsub -A %s -n %d -t %d -q %s --mode script --disable_preboot %s' % (acct, nodes, time, queue, ' '.join(sys.argv))
    print qsub_args
    os.system(qsub_args)
    sys.exit(1)

partsize = int(os.environ['COBALT_PARTSIZE'])
partition = os.environ['COBALT_PARTNAME']
jobid = int(os.environ['COBALT_JOBID'])

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
for job, (block, corner, shape) in zip(jobdirs, block_corner_iter(blocks, npj)):
    print job, (block, corner, shape)
    os.chdir(os.path.join(scratch, job))
    log = open('%d.vasp.stdout' % jobid, 'w')

    locargs = '--block %s --corner %s --shape %s' % (block, corner, shape)
    runjob_args = ('runjob %s -n %d -p %d %s : %s' % (locargs, npj*ppn, ppn, envargs, vasp)).split()
    print ' '.join(runjob_args)
    print

    jobs.append(subprocess.Popen(runjob_args, stdout=log))
    logs.append(log)
    
# wait for all background jobs to finish, then flush their logs
for (job, log) in zip(jobs, logs):
    job.wait()
    log.flush()

