#!/usr/bin/env python
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

import os
import glob
import sys
import time

sys.path.insert(0, os.getcwd())

import numpy as np

from ase.io import read
from ase.constraints import FixAtoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import ase.units as units

from quippy.atoms  import Atoms
from quippy.io     import AtomsWriter
from quippy.farray import fzeros, frange, unravel_index, farray
from quippy import set_fortran_indexing, calc_nye_tensor 
from quippy.system import verbosity_push, PRINT_VERBOSE, enable_timing, system_timer

from quippy.potential import ForceMixingPotential, Potential
from quippy.lotf import LOTFDynamics, update_hysteretic_qm_region

sys.path+=['.']
import params

set_fortran_indexing=False

try:
  if not params.classical:
    from atomsserver import QUIPClient, VaspClient
    from distribfm import DistributedForceMixingPotential
except:
  print 'No BGQ'

def log_pred_corr_errors(dynamics, logfile):
    logline = '%s err %10.1f%12.6f%12.6f\n' % (dynamics.state_label,
                                               dynamics.get_time()/units.fs,
                                               dynamics.rms_force_error,
                                               dynamics.max_force_error)
    print logline
    logfile.write(logline)

def set_quantum(x, n_quantum):
    # parameters for the simulation
    qr = 1       # number of quantum regions per zone
    mom = [3.0 for at in range(len(x))]
    x.set_initial_magnetic_moments(mom)
    # add properties
    hyb = np.zeros(x.n)
    #x.add_property('hybrid_vec', hyb, overwrite=True)
    x.add_property('hybrid', 0, overwrite=True)
    x.add_property('hybrid_vec', 0, overwrite=True)
    x.add_property('hybrid_1', 0)
    x.add_property('hybrid_mark_1', 0)
    core = x.params['core']
    for j in frange(1, x.n):
        if (x.diff_min_image(core,j)[0]**2 + x.diff_min_image(core,j)[1]**2)**0.5 < 5.0 :
            x.hybrid_vec[j-1] = 1
            x.hybrid[j-1] = 1

    x.hybrid_1[:] = x.hybrid_vec[:]  
    x.params['core'] = core[:]
    return x

# Print some information every time step
def printstatus():
    if dynamics.nsteps == 1:
        print """
State      Time/fs    Temp/K
--------------------------------"""

    log_format = ('%(label)-4s%(time)12.1f%(temperature)12.6f')

    atoms.info['label'] = dynamics.state_label  # Label for the status line
    atoms.info['time'] = dynamics.get_time()/units.fs
    atoms.info['temperature'] = (atoms.get_kinetic_energy() /
                                 (1.5*units.kB*len(atoms)))
    print log_format % atoms.info

def traj_writer(dynamics):
    if params.extrapolate_steps == 1 or dynamics.state == LOTFDynamics.Interpolation:
        trajectory.write(dynamics.atoms)

def update_qm_region(atoms):
    cut =3.
    rr = 10
    core = fzeros(3)
    qr = 1 
    thr = 0.1
    core[:] = atoms.params['core']
    print core[1], core[2], 'CORE'
    fixed_mask = (np.sqrt((atoms.positions[:,0]-core[1])**2 + (atoms.positions[:,1]-core[2])**2) < rr)
    cl = atoms.select(mask=fixed_mask, orig_index=True) 
    print 'Number of Atoms in Cluster', cl.n
    cl.set_cutoff(cut)
    cl.calc_connect()
    cl = Atoms(cl)
    x0 = Atoms('ref_slab.xyz')
    x0.set_cutoff(cut)
    x0.calc_connect()
    alpha = calc_nye_tensor(cl, x0, 3, 3, cl.n)    
    cl.screw = alpha[2,2,:]
    cl.edge = alpha[2,0,:]
    sum = 0
    c=farray([0.,0.,0.])
    for i in range(cl.n):
        sum = sum + cl.edge[i]
        c[1] = c[1] + cl.pos[1,i]*cl.edge[i]
        c[2] = c[2] + cl.pos[2,i]*cl.edge[i]
    c[1] = c[1]/sum
    c[2] = c[2]/sum
    c[3] = atoms.lattice[2,2]/2.
    core[:] = c.copy()
    core  = np.array([98.0, 98.0, 1.49]) 
    print 'New Core', core
    old_qm_list = atoms.hybrid_vec.nonzero()[0]
    new_qm_list = update_hysteretic_qm_region(atoms, old_qm_list, core[:],
                                              params.qm_inner_radius,
                                              params.qm_outer_radius,
                                              update_marks=False)
#Force Mixing Potential requires hybrid property:
    atoms.hybrid[:] = 0
    atoms.hybrid[new_qm_list] = 1
#Distributed Force Mixing Properties:
    atoms.hybrid_vec[:] = 0
    atoms.hybrid_vec[new_qm_list] = 1
    atoms.hybrid_1[:] = atoms.hybrid_vec[:]
    atoms.params['core'] = core[:]
    return 

if __name__=='__main__':
  if hasattr(params, 'do_verbose') and params.do_verbose:
    verbosity_push(PRINT_VERBOSE)
  if hasattr(params, 'do_timing') and params.do_timing:
    enable_timing()
# ********** Read input file ************
  input_file = params.input_file
  print 'Loading atoms from file %s' % input_file
  atoms = read(input_file)
  
  if params.continuation:
      # restart from last frame of most recent trajectory file
      traj_files = sorted(glob.glob('[0-9]*.traj.xyz'))
      if len(traj_files) > 0:
          last_traj = traj_files[-1]
          input_file = last_traj + '@-1'
  
  # loading reference configuration for Nye tensor evaluation
  x0 = read(params.reference_file)
  
  # convert to quippy Atoms - FIXME in long term, this should not be necesary
  atoms = Atoms(atoms)
  x0 = Atoms(x0)
  
  x0.set_cutoff(3.0)
  x0.calc_connect()
  
  # ***** Setup constraints *******
  
  
  # ******* Set up potentials and calculators ********
  
  system_timer('init_fm_pot')
  
  mm_pot = Potential(params.mm_init_args,
                     param_filename=params.param_file,
                     cutoff_skin=params.cutoff_skin)
  
  cluster_args = params.cluster_args.copy()
  
  if params.test_mode:
      # dummy QM potential made by swapping Ni and Al species in MM potential               
      qm_pot = Potential(params.mm_init_args,
                         param_filename='m2004flipNiAl.xml',
                         cutoff_skin=params.cutoff_skin)
      qm_clients = qm_pot
      
      # convergence of EAM forces with buffer size is suprisingly slow,
      # so we need to use a big buffer to avoid messing up predictor/corrector
      # error plot
      cluster_args['hysteretic_buffer_inner_radius'] = 12.0
      cluster_args['hysteretic_buffer_outer_radius'] = 14.0
  else:
      qm_clients = params.qm_clients
  
  if params.extrapolate_steps == 1:
      force_mixing_method = 'conserve_momentum'
  else:
      force_mixing_method = 'lotf_adj_pot_minim'
  
  if params.classical:
      # Use the classical EAM only
      atoms.set_calculator(mm_pot)
  else:
      # Use the force mixing potential as the Atoms' calculator
      qmmm_pot = DistributedForceMixingPotential(mm_pot, qm_clients,
                                                 ip=params.ip, port=0,
                                                 rundir=params.rundir,
                                                 cutoff_skin=params.cutoff_skin,
                                                 cluster_args=cluster_args,
                                                 method=force_mixing_method,
                                                 fit_hops=2,
                                                 lotf_spring_hops=2,
                                                 test_mode=params.test_mode,
                                                 save_clusters=params.save_clusters,
                                                 force_restart=params.force_restart)
                                                 
      atoms.set_calculator(qmmm_pot)
  system_timer('init_fm_pot')
  
  if not params.continuation and not params.classical:
      print 'Finding initial dislocation core positions...'
      atoms  = set_quantum(atoms, params.n_core)
      print 'done.'
  
  # ********* Setup and run MD ***********
  # Set the initial temperature to 2*simT: it will then equilibriate to
  # simT, by the virial theorem
  if params.test_mode:
      np.random.seed(0) # use same random seed each time to be deterministic 
  if params.rescale_velo:
      MaxwellBoltzmannDistribution(atoms, 2.0*params.sim_T)
  
  # Save frames to the trajectory every `traj_interval` time steps
  # but only when interpolating
  trajectory = AtomsWriter(os.path.join(params.rundir, params.traj_file))
  # Initialise the dynamical system
  system_timer('init_dynamics')
  check_force_error = params.check_force_error
  if params.extrapolate_steps == 1:
      dynamics = VelocityVerlet(atoms, params.timestep)
      check_force_error = False
      if not params.classical:
          qmmm_pot.set(calc_weights=True)
      dynamics.state_label = 'D'
  else:
      dynamics = LOTFDynamics(atoms, params.timestep,
                              params.extrapolate_steps,
                              check_force_error=check_force_error)
  system_timer('init_dynamics')
  
  # Function to update the QM region at the beginning of each extrapolation cycle   
  if not check_force_error:
      if params.extrapolate_steps == 1:
          if not params.classical:
              dynamics.attach(update_qm_region, 1, dynamics.atoms)
      else:
          dynamics.set_qm_update_func(update_qm_region)
  
  if check_force_error:
      pred_corr_logfile = open(os.path.join(params.rundir, 'pred-corr-error.txt'), 'w')
      dynamics.attach(log_pred_corr_errors, 1, dynamics, pred_corr_logfile)
     
  dynamics.attach(traj_writer, params.traj_interval, dynamics)
  print 'Cutoff of atoms is', dynamics.atoms.cutoff, 'A'
  dynamics.attach(printstatus)
  # Start running!
  system_timer('dynamics_run')
  dynamics.run(params.nsteps)
  system_timer('dynamics_run')
  
  if check_force_error:
      pred_corr_logfile.close()
  
# Shutdown clients and cleanup server
  if not params.classical:
    qmmm_pot.server.shutdown()

