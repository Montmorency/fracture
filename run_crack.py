import os 
import numpy as np
import ase.units as units
from  ase.constraints             import FixAtoms
from  ase.md.verlet               import VelocityVerlet
from  ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from  quippy           import Atoms
from  quippy.potential import Potential
from  quippy.io        import AtomsWriter, AtomsReader
from  quippy.crack     import get_strain, get_energy_release_rate,\
                             ConstantStrainRate, find_crack_tip_stress_field
from quippy           import set_fortran_indexing
from quippy.potential import ForceMixingPotential
from quippy.lotf      import LOTFDynamics, update_hysteretic_qm_region
from simulate_crack   import update_qm_region_context, fix_edges, set_qmmm_pot, pass_print_context,\
                             check_if_cracked_context, pass_trajectory_context
import argparse
#simulation parameters
extrapolate_steps = 10         # Number of steps for predictor-corrector
                               # interpolation and extrapolation
sim_T       = 300.0*units.kB   # Simulation temperature
nsteps      = 5000             # Total number of timesteps to run for
timestep    = 2.0*units.fs     # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang    # Amount by which potential cutoff is increased
                               # for neighbour calculations
tip_move_tol = 12.0            # Distance tip has to move before crack
                               # is taken to be running
#strain_rate    = 0.0*(1.0/units.fs)
strain_rate    = 0.0
traj_interval  = 5             # Number of time steps between interpolations
print_interval = 10             # time steps between trajectory prints 10 fs
param_file     = 'params.xml'   # Filename of XML file containing
                                # potential parameters
mm_init_args = 'IP EAM_ErcolAd' # Classical potential
qm_init_args = 'TB DFTB'        # Initialisation arguments for QM potential
input_file   = 'crack.xyz'      # crack_slab
traj_file    = 'crack_traj.xyz' # Trajectory output file in (NetCDF format)
pot_dir      = '/users/k1511981/pymodules/imeall/imeall/potentials'
#pot_file     = os.path.join(pot_dir, 'Fe_Mendelev.xml')
pot_file     = os.path.join(pot_dir, 'PotBH.xml')

# Restart from last point in trajectory file:
if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input_file",  help="File containing initial configuration. \
                                                   Default is crack.xyz.", default='crack.xyz')
  parser.add_argument("-o", "--output_file", help="File trajectory is written to.", default='crack_traj.xyz')
  parser.add_argument("-p", "--pot_file",    help="Potential is set.", default=os.path.join(pot_dir,'PotBH.xml'))
  parser.add_argument("-r", "--restart",     help="If false thermalizes atoms and fixes boundary conditions,\
                                                   frame in trajectory.", default=False)
  args = parser.parse_args()
  print args.input_file
  restart = args.restart
  print restart
  if restart:
    print 'Restarting job'
  input_file = args.input_file
  pot_file   = args.pot_file
  if args.output_file:
    traj_file = args.output_file
  print 'Output_file: ', traj_file
  print 'POT FILE', pot_file
#Need to add an arg parser here
  learn_on_the_fly = False
  if learn_on_the_fly:
    print 'Initialize Potentials'
    mm_pot = Potential(mm_init_args, param_filename=pot_file, cutoff_skin=cutoff_skin)
    qm_pot = Potential(qm_init_args, param_filename='params.xml')
    print 'Read last MD snapshot'
    print 'Reading from ', input_file
    atoms = Atoms(input_file)
    strain_atoms = fix_edges(atoms)
    current_crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot)
    qmmm_pot = set_qmmm_pot(atoms, current_crack_pos)
    print 'Setup dynamics'
#If input_file is crack.xyz the cell has not been thermalized yet.
#Otherwise it will recover temperature from the previous run.
    if not restart: 
      print 'Thermalizing atoms'
      #Use same random seed so that simulations are deterministic
      np.random.seed(42)
      MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)
    print 'Initializing Dynamics'
    dynamics = LOTFDynamics(atoms, timestep, extrapolate_steps)
    print 'Attaching check_if_cracked dynamics'
    dynamics.attach(check_if_cracked_context(strain_atoms), 1, atoms)
    dynamics.set_qm_update_func(update_qm_region_context(qmmm_pot, atoms))
    print 'Attaching trajectories to dynamics'
    trajectory = AtomsWriter(traj_file)
#Only wriates trajectory if the system is in the LOTFDynamicas
#Interpolation 
    dynamics.attach(pass_trajectory_context(trajectory, dynamics), traj_interval, dynamics)
    dynamics.attach(pass_print_context(atoms, dynamics))
    print 'Running Crack Simulation'
    dynamics.run(nsteps)
    print 'Crack Simulation Finished'
  else:
    mm_pot = Potential(mm_init_args, param_filename=pot_file, cutoff_skin=cutoff_skin)
    atoms = AtomsReader(input_file)[-1]
    atoms = Atoms(atoms)
    atoms.set_calculator(mm_pot)
    strain_atoms = fix_edges(atoms)
    current_crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot)
    print 'Current Crack Position: ', current_crack_pos
    if not restart: 
      print 'Thermalizing Atoms'
      MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)
    dynamics = VelocityVerlet(atoms, timestep)
    dynamics.attach(pass_print_context(atoms, dynamics, mm_pot))
    dynamics.attach(check_if_cracked_context(strain_atoms, mm_pot), 1, atoms)
    trajectory = AtomsWriter(traj_file)
    dynamics.attach(pass_trajectory_context(trajectory, dynamics), print_interval, dynamics)
    print 'Running Crack Simulation'
    dynamics.run(nsteps)
    print 'Crack Simulation Finished'
