#!/usr/bin/env python
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)

import os
import glob
import sys
import time

import argparse

sys.path.insert(0, os.getcwd())

import numpy as np

from ase.constraints import FixAtoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import ase.units as units

from quippy           import set_fortran_indexing, calc_nye_tensor 
from quippy.io        import AtomsWriter
from quippy.lotf      import LOTFDynamics, update_hysteretic_qm_region
from quippy.crack     import get_strain, get_energy_release_rate,\
                             ConstantStrainRate, find_crack_tip_stress_field
from quippy.atoms     import Atoms
from quippy.system    import verbosity_push, PRINT_VERBOSE, enable_timing, system_timer
from quippy.potential import ForceMixingPotential, Potential
from quippy.clusters  import  HYBRID_NO_MARK, HYBRID_ACTIVE_MARK
from simulate_crack   import update_qm_region_context, fix_edges, set_qmmm_pot, pass_print_context,\
                             check_if_cracked_context, pass_trajectory_context

from matscipy.socketcalc  import  VaspClient, SocketCalculator
from bgqtools             import (get_bootable_blocks, boot_blocks, block_corner_iter, 
                                  get_hostname_ip, get_cobalt_info, set_unbuffered_stdout)
from distribfm            import DistributedForceMixingPotential

sys.path+=['.']

set_fortran_indexing(False)


def log_pred_corr_errors(dynamics, logfile):
    logline = '%s err %10.1f%12.6f%12.6f\n' % (dynamics.state_label,
                                               dynamics.get_time()/units.fs,
                                               dynamics.rms_force_error,
                                               dynamics.max_force_error)
    print logline
    logfile.write(logline)

def set_quantum_disloc(x):
    # parameters for the simulation
    mom = [3.0 for at in range(len(x))]
    x.set_initial_magnetic_moments(mom)
    # add properties
    x.add_property('hybrid', 0, overwrite=True)
    x.add_property('hybrid_vec', 0, overwrite=True)
    x.add_property('hybrid_1', 0)
    x.add_property('hybrid_mark_1', 0)
    core = x.params['core']
    for j in range(len(x)):
    #Function calls still require fortran indexing:
        if (x.diff_min_image(core,j)[1]**2 + x.diff_min_image(core,j)[2]**2)**0.5 < 5.0 :
            x.hybrid_vec[j] = 1
            x.hybrid[j]     = 1
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
    if extrapolate_steps == 1 or dynamics.state == LOTFDynamics.Interpolation:
        trajectory.write(dynamics.atoms)

def update_qm_region(atoms, dis_type='edge'):
    cut =3.
    rr = 10
    qr = 1 
    thr = 0.1
    core[:] = atoms.params['core']
    fixed_mask = (np.sqrt((atoms.positions[:,0]-core[0])**2 + (atoms.positions[:,1]-core[1])**2) < rr)
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
    cl.edge  = alpha[2,0,:]
    if dis_type  == 'screw':
      defect_pos = cl.screw
    elif dis_type == 'edge':
      defect_pos = cl.edge
    total_def = 0.0
    c = np.array([0.,0.,0.])
    mom = [3.0 for at in range(len(atoms))]
    atoms.set_initial_magnetic_moments(mom)
    for i in range(cl.n):
        defect_pos = defect_pos   + cl.edge[i]
        c[1] = c[1] + cl.positions[i,0]*defect_pos[i]
        c[2] = c[2] + cl.pos[i,0]*defect_pos[i]
    c[0] = c[0]/total_def
    c[1] = c[1]/total_def
    c[2] = atoms.lattice[2,2]/2.
    core[:] = c.copy()
    old_qm_list = atoms.hybrid_vec.nonzero()[0]
    new_qm_list = update_hysteretic_qm_region(atoms, old_qm_list, core[:],
                                              qm_inner_radius,
                                              qm_outer_radius,
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

def update_qm_region_crack(atoms):
  mm_pot = Potential(mm_init_args,
                     param_filename = param_file,
                     cutoff_skin    = cutoff_skin)
  crack_pos    = find_crack_tip_stress_field(atoms, calc=mm_pot)
  old_qm_list  = atoms.hybrid_vec.nonzero()[0]
  new_qm_list  = update_hysteretic_qm_region(atoms, old_qm_list, crack_pos, 
                                             qm_inner_radius, 
                                             qm_outer_radius, 
                                             update_marks=False)
#lets try setting the atoms object properties by hand?
  atoms.hybrid[:] = 0
  atoms.hybrid[new_qm_list] = 1
#Distributed Force Mixing Properties:
  atoms.hybrid_vec[:] = 0
  atoms.hybrid_vec[new_qm_list] = 1
  atoms.hybrid_1[:] = atoms.hybrid_vec[:]
  atoms.params['core'] = crack_pos
  atoms.params['CrackPos'] = crack_pos
  return

if __name__=='__main__':
  if do_verbose:
    verbosity_push(PRINT_VERBOSE)
  if do_timing:
    enable_timing()
# ********** Simulation Arguments  ************ #
  parser = argparse.ArgumentParser()
  parser.add_argument("-g","--geom", default='crack')
  parser.add_argument("-inp", "--input_file", required=True)
  parser.add_argument("-st",  "--sim_T", help='Simulation Temperature in Kelvin. Default is 300 K.', type=float, required=True)
  parser.add_argument("-cfe", "--check_force_error", help='Perform a DFT calculation at each step in the trajectory.', action='store_true')
# ********** Parallel Arguments  ************ #
  parser.add_option("--block")
  parser.add_option("--corner")
  parser.add_option("--shape")
  parser.add_option("--npj")
  parser.add_option("--ppn")

# parse args string:
  args        = parser.parse_args()

  if args.input_file == '':
    input_file = 'crack.xyz'
  else:
    input_file = args.input_file

  sim_T             = args.sim_T*units.kB
  geom              = args.geom
  check_force_error = args.check_force_error
  continuation      = False

# COBALT CONFIGURATION
  vasp  = '/projects/SiO2_Fracture/iron/vasp.bgq'
  block  = args.block
  corner = args.corner
  shape  = args.shape
  npj    = int(args.npj)
  ppn    = int(args.ppn)
  hostname, ip = get_hostname_ip()

  do_timing    = False
  do_verbosity = False

  reference_file = 'ref_slab.xyz' # Reference file for Nye tensor
  continuation = False             # If true, restart form last frame of most recent *.traj.xyz file
  classical    = False             # If true, do classical MD instead of QM/MM
  sim_T        = 300.0*units.kB    # Simulation temperature
  rescale_velo = False             # Rescale velocities to 2*sim_T  
  timestep     = 1.0*units.fs      # Timestep (NB: time base units are not fs!)
  cutoff_skin  = 2.0*units.Ang     # Amount by which potential cutoff is increased
                                 # for neighbour calculations
  traj_file    = '%s.traj.xyz'        # Trajectory output file

  try:
    pot_dir      = os.environ['POTDIR']
  except:
    sys.exit("POTDIR variable not set in bash environment.")

  param_file   = 'PotBH.xml'   # Filename of XML file containing
  param_file   = os.path.join(pot_dir, param_file)
                                   # potential parameters
  mm_init_args = 'IP EAM_ErcolAd'  # Initialisation arguments for
                                   # classical potential
# ADDITIONAL PARAMETERS FOR THE QM/MM simulation:
  qm_inner_radius   = 3.0*units.Ang # Inner hysteretic radius for QM region
  qm_outer_radius   = 5.0*units.Ang # Outer hysteretic radius for QM region
  hyst_buffer_inner = 7.0 # Inner hysteretic radius for QM region
  hyst_buffer_outer = 9.0 # Outer hysteretic radius for QM region
  extrapolate_steps = 5   # Number of steps for predictor-corrector
                          # interpolation and extrapolation
  traj_interval     = extrapolate_steps # Only print the DFT steps
  check_force_error = False         # Do QM calc at each step to check pred/corr. err
  nsteps = 1000
  
  save_clusters = False  # if True, write each QM cluster to .xyz files
  force_restart = True   # if True, force VASP to restart for each QM cluster
  
  traj_index = 1
  traj_file = '%d.traj.xyz' % traj_index
  while os.path.exists(traj_file):
    traj_index += 1
    traj_file = '%d.traj.xyz' % traj_index

  cluster_args = dict(single_cluster=False,
                      cluster_calc_connect=False,
                      cluster_hopping=False,
                      cluster_hopping_nneighb_only=True,
                      cluster_periodic_z = True, # match Gamma vs. kpts
                      cluster_vacuum     = 7.0,
                      hysteretic_buffer=True,
                      hysteretic_buffer_inner_radius=hyst_buffer_inner,
                      hysteretic_buffer_outer_radius=hyst_buffer_outer,
                      min_images_only=True,
                      terminate=False,
                      force_no_fix_termination_clash=True,
                      randomise_buffer=False)

  vasp_args=dict(parmode='cobalt', xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001, 
                 kpts=[1, 1, 8], kpar=4, lreal='auto', ibrion=13, nsw=1000000, nelmdl=-15, ispin=2,
                 nelm=100, algo='VeryFast', npar=32, lplane=False, lwave=False, lcharg=False, istart=0,
                 voskown=1, ismear=1, sigma=0.1, isym=0) # possibly try iwavpr=12, should be faster if it works

  para_args  = dict(npj=npj, ppn=ppn, block=block, corner=corner, shape=shape, exe=vasp)
  all_args   = dict(vasp_args.items() + para_args.items())
  qm_clients = [VaspClient(ip=ip, bgq=True, **all_args)] 

##### Finished VASP INITIALIZATION AND SOCKET CONFIGURATION ######
##### MAIN PROGRAM ############

  print 'Loading atoms from file %s' % input_file
  atoms = Atoms(input_file)
  atoms = Atoms(atoms)
  if continuation:
      # restart from last frame of most recent trajectory file
      traj_files = sorted(glob.glob('[0-9]*.traj.xyz'))
      if len(traj_files) > 0:
          last_traj = traj_files[-1]
          input_file = last_traj + '@-1'
  
  # loading reference configuration for Nye tensor evaluation
  # convert to quippy Atoms - FIXME in long term, this should not be necesary
  x0 = Atoms(params.reference_file)
  x0 = Atoms(x0)
  x0.set_cutoff(3.0)
  x0.calc_connect()
  
  print params.param_file
  # ******* Set up potentials and calculators ********
  system_timer('init_fm_pot')
  pot_file = os.path.join(params.pot_dir, params.param_file)
  mm_pot = Potential(mm_init_args,
                     param_filename=params.param_file,
                     cutoff_skin=params.cutoff_skin)
  
  cluster_args = params.cluster_args.copy()
  
  if test_mode:
      # dummy QM potential made by swapping Ni and Al species in MM potential               
      qm_pot = Potential(mm_init_args,
                         param_filename='m2004flipNiAl.xml',
                         cutoff_skin=params.cutoff_skin)
      qm_clients = qm_pot
      
      # convergence of EAM forces with buffer size is suprisingly slow,
      # so we need to use a big buffer to avoid messing up predictor/corrector
      # error plot
      cluster_args['hysteretic_buffer_inner_radius'] = 12.0
      cluster_args['hysteretic_buffer_outer_radius'] = 14.0
  else:
  #    qm_clients = params.qm_clients
    pass
  
  if params.extrapolate_steps == 1:
      force_mixing_method = 'conserve_momentum'
  else:
      force_mixing_method = 'lotf_adj_pot_minim'
  
      # Use the classical EAM only
      atoms.set_calculator(mm_pot)
  else:
      # Use the force mixing potential as the Atoms' calculator
      qmmm_pot = DistributedForceMixingPotential(mm_pot, qm_clients,
                                                 ip               = ip, 
                                                 port             = 0,
                                                 rundir           = rundir,
                                                 cutoff_skin      = cutoff_skin,
                                                 cluster_args     = cluster_args,
                                                 method           = force_mixing_method,
                                                 fit_hops         = 2,
                                                 lotf_spring_hops = 2,
                                                 test_mode        = test_mode,
                                                 save_clusters    = save_clusters,
                                                 force_restart    = force_restart)
                                                 
      atoms.set_calculator(qmmm_pot)
  system_timer('init_fm_pot')
  
  if not continuation and not classical:
      print 'Finding initial Quantum Core positions...'
      if geom =='disloc':
        atoms  = set_quantum_disloc(atoms)
      elif geom=='crack':
        mom   = [3.0 for at in range(len(atoms))]
        atoms.set_initial_magnetic_moments(mom)
        atoms.add_property('hybrid', 0, overwrite=True)
        atoms.add_property('hybrid_vec', 0, overwrite=True)
        atoms.add_property('hybrid_1', 0)
        atoms.add_property('hybrid_mark_1', 0)
        crackpos = atoms.info['CrackPos']
        qm_list_old = []
        qm_list   = update_hysteretic_qm_region(atoms, [], crackpos, qm_inner_radius,
                                                qm_outer_radius,
                                                update_marks=False)
        atoms.hybrid[qm_list]     = HYBRID_ACTIVE_MARK
        atoms.hybrid_vec[qm_list] = HYBRID_ACTIVE_MARK
        atoms.hybrid_1[qm_list]   = HYBRID_ACTIVE_MARK
        atoms.hybrid_mark_1[qm_list] = HYBRID_ACTIVE_MARK
        atoms.params['core'] = crackpos
        print HYBRID_ACTIVE_MARK
        print 'Core Found. No. Quantum Atoms:', sum(atoms.hybrid[:])
      else:
        print 'No cell geometry given, specifiy either disloc or crack', 1/0 
  
  # ********* Setup and run MD ***********
  # Set the initial temperature to 2*simT: it will then equilibriate to
  # simT, by the virial theorem
  if test_mode:
      np.random.seed(0) # use same random seed each time to be deterministic 
  if rescale_velo:
      MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)
  # Save frames to the trajectory every `traj_interval` time steps
  # but only when interpolating
  trajectory = AtomsWriter(os.path.join(rundir, traj_file))
  # Initialise the dynamical system
  system_timer('init_dynamics')
  if extrapolate_steps == 1:
      dynamics = VelocityVerlet(atoms, timestep)
      check_force_error = False
      if not classical:
          qmmm_pot.set(calc_weights=True)
      dynamics.state_label = 'D'
  else:
    print 'Initializing LOTF Dynamics'
    dynamics = LOTFDynamics(atoms, timestep,
                            extrapolate_steps,
                            check_force_error=check_force_error)
  system_timer('init_dynamics')
#Function to update the QM region at the beginning of each extrapolation cycle   
  if not check_force_error:
      if extrapolate_steps == 1:
          if not classical:
              dynamics.attach(update_qm_region, 1, dynamics.atoms)
      else:
      # choose appropriate update function for defects or crack.
      # or grainboundary.
        print 'Setting Update Function'
        if geom =='disloc':
          dynamics.set_qm_update_func(update_qm_region)
        elif geom =='crack':
          dynamics.set_qm_update_func(update_qm_region_crack)
        else:
          print 'No geometry chosen', 1/0
  
  if check_force_error:
      pred_corr_logfile = open(os.path.join(rundir, 'pred-corr-error.txt'), 'w')
      dynamics.attach(log_pred_corr_errors, 1, dynamics, pred_corr_logfile)
     
  dynamics.attach(traj_writer, traj_interval, dynamics)
  print 'Cutoff of atoms is ', dynamics.atoms.cutoff, 'A'
  dynamics.attach(printstatus)
  # Start running!
  system_timer('dynamics_run')
  dynamics.run(nsteps)
  system_timer('dynamics_run')
  
  if check_force_error:
      pred_corr_logfile.close()
  
# Shutdown clients and cleanup server
  if not classical:
    qmmm_pot.server.shutdown()

