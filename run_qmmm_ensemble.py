#!/usr/bin/env python
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

import os
import sys
import glob
import socket
import shutil
import argparse
import subprocess

from bgqtools import (get_bootable_blocks, boot_blocks, block_corner_iter,
                      set_unbuffered_stdout)

set_unbuffered_stdout()


parser = argparse.ArgumentParser()

parser.add_argument("-t",   "--time", type=int, default = 60)
parser.add_argument("-jp",  "--jobpattern",     default = 'T3*')
parser.add_argument("-js",  "--jobstart", type=int, default  = 0)
parser.add_argument("-jf",  "--jobfinish", type=int, default = 8)
parser.add_argument("-npj", "--npj",  type=int, default = 128)
parser.add_argument("-ppn", "--ppn",  type=int, default =  2)
parser.add_argument("-c",   "--continuation", action='store_true')

args       = parser.parse_args()
jobstart   = args.jobstart
jobfinish  = args.jobfinish
jobpattern = args.jobpattern
npj        = args.npj # nodes per job
ppn        = args.ppn # MPI tasks per node

if args.continuation:
  cont_string = '-c'
else:
  cont_string = ''
  

acct  = 'SiO2_Fracture'
time  = args.time
queue = 'default'
mapping = 'ABCDET'
scratch = os.getcwd()
vasp    = '/projects/SiO2_Fracture/iron/vasp.bgq'
envargs = '--envs RUNJOB_MAPPING=%s --envs MPIRUN_ENABLE_TTY_REPORTING=0' % mapping

hostname = socket.gethostname()
print 'Hostname: %s' % hostname

#The submit resubmit pattern means the shell will strip the quotes from the
#job pattern in order to not have wild cards expanded in the shell we need the quotes
#by double enclosing the quotes we protect it from the wildcard expand then
#we strip the extra quotes to find the directors and determine the nodes
#then when subprocess resubmits this script the quotes are again included so when qsub
#executes the script we have the pattern enclosed in quotes '"pattern"'
jobdirs = glob.glob(jobpattern.replace("\"",''))[jobstart:jobfinish]

print 'jobdirs = %s' % jobdirs

njobs = len(jobdirs)
nodes = npj*njobs

if 'COBALT_PARTSIZE' not in os.environ:
    print 'Not running under control of cobalt. Launching qsub...'
    qsub_args = 'qsub -A {acct} -n {nodes} -t {time} -q {queue} --mode script --disable_preboot {script} {flags}'.format(acct=acct, 
                        nodes=nodes, time=time, queue=queue, script=sys.argv[0], flags=' '.join(sys.argv[1:]))
    print qsub_args
    job = subprocess.Popen(qsub_args.split())
    job.wait()
    sys.exit(1)

partsize  = int(os.environ['COBALT_PARTSIZE'])
partition = os.environ['COBALT_PARTNAME']
jobid     = int(os.environ['COBALT_JOBID'])

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
    shutil.copy(os.path.join('/home/lambert/pymodules/fracture', 'run_qmmm_e.py'), os.path.join(scratch,job))
    os.chdir(os.path.join(scratch, job))
    log = open('%d.vasp.stdout' % jobid, 'w')
# We pass lots of arguments to the run_qmmm script that tell it what to do
# sim temperature, check force error etc.
# and where (i.e. which block corner shape) to do it.
    physargs    = '-st {0} {1}'.format(300, cont_string)
    locargs     = '--block {0} --corner {1} --shape {2}'.format(block, corner, shape)
    nodeargs    = '--npj {0} --ppn {1}'.format(npj, ppn)
    runjob_args = ('python run_qmmm_e.py {physargs} {locargs} {nodeargs} '.format(physargs=physargs, locargs=locargs, nodeargs=nodeargs)).split()
    jobs.append(subprocess.Popen(runjob_args, stdout=log))
    logs.append(log)
# wait for all background jobs to finish, then flush their logs
for (job, log) in zip(jobs, logs):
    job.wait()
    log.flush()
