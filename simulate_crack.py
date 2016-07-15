import numpy as np
from   ase.constraints import FixAtoms
from   ase.md.verlet import VelocityVerlet
from   ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import ase.units as units

from quippy import Atoms
from quippy.potential import Potential
from quippy.io import AtomsWriter, AtomsReader
from quippy.crack import(get_strain, get_energy_release_rate,
                         ConstantStrainRate,
                         find_crack_tip_stress_field)
from quippy import set_fortran_indexing
from quippy.potential import ForceMixingPotential
from quippy.lotf import LOTFDynamics, update_hysteretic_qm_region
from quippy.clusters import HYBRID_NO_MARK, HYBRID_ACTIVE_MARK

#simulation parameters
qm_init_args      = 'TB DFTB'       # Initialisation arguments for QM potential
qm_inner_radius   = 18.0*units.Ang   # Inner hysteretic radius for QM region
qm_outer_radius   = 21.0*units.Ang  # Inner hysteretic radius for QM region
extrapolate_steps = 10         # Number of steps for predictor-corrector
                               # interpolation and extrapolation
input_file  = 'crack.xyz'      # crack_slab
sim_T       = 300.0*units.kB   # Simulation temperature
nsteps      = 10000            # Total number of timesteps to run for
timestep    = 1.0*units.fs     # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang    # Amount by which potential cutoff is increased
                               # for neighbour calculations
tip_move_tol = 10.0            # Distance tip has to move before crack
                               # is taken to be running
#Need have simulate crack use the same crack dictionary so I don't need to do this by hand.
#strain_rate = 1e-5*(1.0/units.fs) 
strain_rate    = 1.e-5*(1.0/(7.93*units.fs))
traj_interval  = 10             # Number of time steps between
traj_file      = 'traj_lotf_2.xyz'    # Trajectory output file in (NetCDF format)
restart_traj_file      = 'traj_lotf_2b.xyz'    # Trajectory output file in (NetCDF format)
print_interval = 10            # time steps between trajectory prints 10 fs
param_file = 'params.xml'      # Filename of XML file containing
                               # potential parameters
mm_init_args = 'IP SW'         # Classical potential
set_fortran_indexing(False)

#sw_pot = Potential('IP SW', param_filename='params.xml', cutoff_skin=cutoff_skin)
#qm_pot = Potential(qm_init_args, param_filename='params.xml')
restart = True
#if from scratch we load the original cell:
#atoms = AtomsReader('crack.xyz')
#Some general procedures we use to setup crack cells (probably duplicating some stuff in crack.py):


def fix_edges(atoms):
  orig_height    = atoms.info['OrigHeight']
  top    = atoms.positions[:,1].max()
  bottom = atoms.positions[:,1].min()
  left   = atoms.positions[:,0].min()
  right  = atoms.positions[:,0].max()
  fixed_mask = ((abs(atoms.positions[:, 1] - top) < 2.0) |
               (abs(atoms.positions[:, 1] - bottom) < 2.0))
  fix_atoms  = FixAtoms(mask=fixed_mask)
  print('Fixed %d atoms\n' % fixed_mask.sum()) # Print the number of fixed atoms
  strain_atoms = ConstantStrainRate(orig_height, strain_rate*timestep)
  atoms.set_constraint([fix_atoms, strain_atoms])
  return strain_atoms

def pass_print_context(atoms, dynamics):
  def printstatus():
    #if dynamics.nsteps == 1:
    if (dynamics.nsteps%10)==0:
      print """
State      Time/fs    Temp/K     Strain      G/(J/m^2)  CrackPos/A D(CrackPos)/A 
---------------------------------------------------------------------------------"""
      log_format = ('%(label)-4s%(time)12.1f%(temperature)12.6f'+
                    '%(strain)12.5f%(G)12.4f%(crack_pos_x)12.2f    (%(d_crack_pos_x)+5.2f)')
      log_format2 = ('%(label)-4s%(time)12.1f%(temperature)12.6f')
      try:
        atoms.info['label'] = dynamics.state_label                # Label for the status line
      except AttributeError:
        atoms.info['label'] = 'classical'                # Label for the status line
        

      atoms.info['time']  = dynamics.get_time()/units.fs
      atoms.info['temperature'] = (atoms.get_kinetic_energy() /
                                   (1.5*units.kB*len(atoms)))
      atoms.info['strain'] = get_strain(atoms)
      atoms.info['G']      = get_energy_release_rate(atoms)/(units.J/units.m**2)
      try:
        orig_crack_pos = atoms.info['CrackPos'].copy()
        crack_pos = find_crack_tip_stress_field(atoms, calc=sw_pot)
        atoms.info['crack_pos_x']   = crack_pos[0]
        atoms.info['d_crack_pos_x'] = crack_pos[0] - orig_crack_pos[0]
        print log_format % atoms.info
      except KeyError:
        print log_format2 % atoms.info
  return printstatus

#pot1 is the low precision
#pot2 is the high precision, i.e. QM potential
def set_qmmm_pot(atoms, crack_pos):
  qmmm_pot = ForceMixingPotential(pot1=sw_pot, pot2=qm_pot, atoms=atoms,
                               qm_args_str='single_cluster cluster_periodic_z carve_cluster '+
                              'terminate cluster_hopping=F randomise_buffer=F',
                               fit_hops=4,
                               lotf_spring_hops=3,
                               hysteretic_buffer=True,
                               hysteretic_buffer_inner_radius=qm_inner_radius,
                               hysteretic_buffer_outer_radius=qm_outer_radius,
                               cluster_hopping_nneighb_only=False,
                               min_images_only=True)
  qmmm_pot.atoms = atoms
  atoms.set_calculator(qmmm_pot)
  qm_list = update_hysteretic_qm_region(atoms, [], crack_pos, qm_inner_radius,
                                        qm_outer_radius, update_marks=True)
  qmmm_pot.set_qm_atoms(qm_list)
  return qmmm_pot

def pass_trajectory_context(trajectory, dynamics):
  def traj_writer(dynamics):
    if dynamics.state == LOTFDynamics.Interpolation:
      trajectory.write(dynamics.atoms)
  return traj_writer

# Prevents Strain from being incremented behind the crack tip
def check_if_cracked_context(strain_atoms):
  def check_if_cracked(atoms):
    orig_crack_pos = atoms.info['CrackPos'].copy()
    crack_pos = find_crack_tip_stress_field(atoms, calc=sw_pot)
    if not atoms.info['is_cracked'] and (crack_pos[0] - orig_crack_pos[0]) > tip_move_tol:
      atoms.info['is_cracked'] = True
      del atoms.constraints[atoms.constraints.index(strain_atoms)]
  return check_if_cracked

def update_qm_region_context(qmmm_pot, atoms):
  def update_qm_region(atoms):
    crack_pos = find_crack_tip_stress_field(atoms, calc=sw_pot)
    qm_list   = qmmm_pot.get_qm_atoms()
    qm_list   = update_hysteretic_qm_region(atoms, qm_list, crack_pos, qm_inner_radius, 
                                            qm_outer_radius, update_marks=True)
    qmmm_pot.set_qm_atoms(qm_list)
#lets try setting the atoms object properties by hand?
    atoms.hybrid[:] = HYBRID_NO_MARK
    atoms.hybrid[qm_list] = HYBRID_ACTIVE_MARK
  return update_qm_region

if __name__=='__main__':
#Randomize initial positions
  MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)
#dynamics = VelocityVerlet(atoms, timestep)
  dynamics = LOTFDynamics(atoms, timestep, extrapolate_steps)
  dynamics.attach(check_if_cracked, 1, atoms)
  dynamics = LOTFDynamics(atoms, timestep, extrapolate_steps)
  dynamics.set_qm_update_func(update_qm_region)
#Writing the trajectory
  trajectory = AtomsWriter(traj_file)
  dynamics.attach(traj_writer, traj_interval, dynamics)
  dynamics.attach(printstatus)
  print 'Running Crack Simulation'
  dynamics.run(nsteps)
  print 'Crack Simulation Finished'
