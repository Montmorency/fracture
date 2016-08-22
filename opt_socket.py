#!/usr/bin/env python
import logging
logging.root.setLevel(logging.DEBUG)

import os
import glob
from distutils    import spawn
from ase.optimize.precon import LBFGS
from   ase.constraints   import UnitCellFilter, StrainFilter, FixedLine

from bgqtools import (get_bootable_blocks, boot_blocks, block_corner_iter,
                      get_hostname_ip, get_cobalt_info,
                      set_unbuffered_stdout)

#from ase import Atoms
from quippy import Atoms, bcc
from matscipy.socketcalc import VaspClient, SocketCalculator

#look for mpirun and vasp on $PATH
#mpirun = spawn.find_executable('mpirun')
#vasp = spawn.find_executable('vasp')
#vasp = '/home/eng/essswb/vasp5/vasp.5.3.new/vasp'
mpirun ='/usr/bin/cobalt-mpirun'
vasp  = '/projects/SiO2_Fracture/iron/vasp.bgq'
npj   = 512 # nodes per job
ppn   = 4
njobs = 1
nodes = npj*njobs

job_xyz = glob.glob('feb*.xyz')[0]
bulk = Atoms(job_xyz)

line=[]
for at in bulk:
  line.append( FixedLine(at.index, (0,0,1)) )
bulk.set_constraint(line)

hostname, ip =  get_hostname_ip()
partsize, partition, job_name = get_cobalt_info()
blocks = get_bootable_blocks(partition, nodes)
print('Available blocks: %s' % blocks)
boot_blocks(blocks)

block, corner, shape = list(block_corner_iter(blocks, npj))[0]
print block, corner, shape

vasp_client = VaspClient(client_id=0,
                         kpts =[16,16,1],
                         amix = 0.01,
                         amin = 0.001,
                         bmix = 0.001,
                         amix_mag = 0.01,
                         bmix_mag = 0.001,
                         npj=npj,
                         ppn=ppn,
                         block=block,
                         corner=corner,
                         shape=shape,
                         exe    = vasp,
                         ediffg = -0.05,
                         nelmin = 6,
                         nelmdl = -15,
                         kpar   = 32,
                         parmode='cobalt',
                         xc='PBE',
                         lreal=False, ibrion=13, nsw=40,
                         algo='VeryFast', npar=8, 
                         lplane=False, lwave=False, lcharg=False, nsim=1,
                         voskown=1, ismear=1, sigma=0.01, iwavpr=11, isym=2, nelm=100)

sock_calc = SocketCalculator(vasp_client, ip=ip, bgq=True)

bulk.set_calculator(sock_calc)
opt = LBFGS(bulk)
opt.run(fmax=0.05, steps=100)

sock_e = bulk.get_potential_energy()
sock_f = bulk.get_forces()
sock_s = bulk.get_stress()

print 'energy', sock_e
print 'forces', sock_f
print 'stress', sock_s

sock_calc.shutdown()
