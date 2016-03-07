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
from quippy.lotf      import LOTFDynamics, update_hysteretic_qm_region

from simulate_crack import pass_qmmm_context, fix_edges, set_qmmm_pot, pass_print_context, \
													 check_if_cracked_context, pass_trajectory_context

#simulation parameters
qm_init_args      = 'TB DFTB'       # Initialisation arguments for QM potential
qm_inner_radius   = 15.0*units.Ang  # Inner hysteretic radius for QM region
qm_outer_radius   = 18.0*units.Ang  # Inner hysteretic radius for QM region
extrapolate_steps = 10         # Number of steps for predictor-corrector
                               # interpolation and extrapolation
sim_T       = 300.0*units.kB   # Simulation temperature
nsteps      = 10000            # Total number of timesteps to run for
timestep    = 1.0*units.fs     # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang    # Amount by which potential cutoff is increased
                               # for neighbour calculations
tip_move_tol = 10.0            # Distance tip has to move before crack
                               # is taken to be running
strain_rate = 1e-5*(1.0/units.fs) 
traj_interval = 10             # Number of time steps between
print_interval = 10            # time steps between trajectory prints 10 fs
param_file = 'params.xml'      # Filename of XML file containing
                               # potential parameters
mm_init_args = 'IP SW'         # Classical potential
input_file  = 'crack.xyz'      # crack_slab
traj_file      = 'traj_lotf_2b.xyz'    # Trajectory output file in (NetCDF format)


#Restart from last point in trajectory file:
if __name__=='__main__':
	print 'Initialize Potentials'
	sw_pot = Potential(mm_init_args, param_filename='params.xml', cutoff_skin=cutoff_skin)
	qm_pot = Potential(qm_init_args, param_filename='params.xml')
	print 'Read last MD snapshot'
	atoms = AtomsReader(input_file)[-1]
	atoms = Atoms(atoms)
	strain_atoms = fix_edges(atoms)
	current_crack_pos = find_crack_tip_stress_field(atoms, calc=sw_pot)

#Not sure how to carry over qmmm atoms from previous calculation:
	qmmm_pot = set_qmmm_pot(atoms, current_crack_pos)
	print 'Setup dynamics'
	dynamics = LOTFDynamics(atoms, timestep, extrapolate_steps)
	dynamics.attach(check_if_cracked_context(strain_atoms), 1, atoms)
	dynamics.set_qm_update_func(pass_qmmm_context(qmmm_pot, atoms))
	print 'Attaching trajectories to dynamics'
	trajectory = AtomsWriter(traj_file)
	dynamics.attach(pass_trajectory_context(trajectory, dynamics), traj_interval, dynamics)
	dynamics.attach(pass_print_context(atoms, dynamics))
	print 'Running Crack Simulation'
	dynamics.run(nsteps)
	print 'Crack Simulation Finished'
