import os 
import argparse
import numpy as np
import ase.units as units

from  ase.constraints             import FixAtoms
from  ase.md.verlet               import VelocityVerlet
from  ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from  ase.optimize     import FIRE, LBFGS
from  ase.lattice.cubic       import BodyCenteredCubic
from  quippy           import Atoms, supercell
from  quippy.potential import Potential
from  quippy.io        import AtomsWriter, AtomsReader
from  quippy.crack     import get_strain, get_energy_release_rate,\
                             ConstantStrainRate, find_crack_tip_stress_field
from quippy           import set_fortran_indexing
from quippy.potential import ForceMixingPotential
from quippy.lotf      import LOTFDynamics, update_hysteretic_qm_region
from fracture.simulate_crack   import update_qm_region_context, fix_edges, set_qmmm_pot, pass_print_context,\
                                      check_if_cracked_context, pass_trajectory_context
from   fracture      import tb_pot
from   quippy.system import verbosity_push, PRINT_VERBOSE


#simulation parameters
set_fortran_indexing(False)
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
try:
  pot_dir      = os.environ['POTDIR']
except:
  sys.exit("PLEASE SET export POTDIR='path/to/potfiles/'")

pot_file     = os.path.join(pot_dir, 'PotBH.xml')

# Restart from last point in trajectory file:
if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input_file",  help="File containing initial configuration. \
                                                   Default is crack.xyz.", default='crackH.xyz')
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
  learn_on_the_fly = True
  if learn_on_the_fly:
    print 'Initialize Potentials'
    atoms = Atoms(input_file)
    cluster = Atoms('cluster.xyz')
    mm_pot = Potential(mm_init_args, param_filename=pot_file, cutoff_skin=cutoff_skin)
    unit_cell = BodyCenteredCubic(directions = [[1,0,0], [0,1,0],[0,0,1]],
                                size = (4,4,4), symbol='Fe', pbc=(1,1,1),
                                latticeconstant = 2.87)

    def proto_qm_pot_callback(unit_cell):
      global qm_pot
      tb = tb_pot.TightBindingPot(alat=1.00, nk=1)
      tb.write_control_file('ctrl.fe', unit_cell)

      qm_pot = Potential('IP LMTO_TBE', param_str="""
<params>
  <LMTO_TBE_params n_types="2" control_file="ctrl.fe">
    <per_type_data type="1" atomic_num="26"/>
    <per_type_data type="2" atomic_num="1"/>
  </LMTO_TBE_params>
</params>""")
      unit_cell.add_property('forces', 0.0, n_cols=3)

#to use tightbinginding
    proto_qm_pot = Potential(callback=proto_qm_pot_callback)
    mm_init_args = 'IP EAM_ErcolAd do_rescale_r=T r_scale=1.00894848312' # Classical potential
    mm_pot = Potential(mm_init_args, param_filename=pot_file, cutoff_skin=cutoff_skin)
    print 'Read last MD snapshot'
    print 'Reading from ', input_file
    TB = False
    if TB:
      qmmm_pot = set_qmmm_pot(atoms, atoms.params['CrackPos'], mm_pot, proto_qm_pot)
# do a one-off calculationto set up the TB potential with correct cluster
      qmmm_pot.get_forces()
    else:
      print 'Running WITH EAM as embedded cluster'
      qm_pot_file  = os.path.join(pot_dir, 'PotBH_fakemod.xml')
      print qm_pot_file
      mm_init_args = 'IP EAM_ErcolAd do_rescale_r=T r_scale=1.01' # Classical potential
      qm_pot       = Potential(mm_init_args, param_filename=qm_pot_file, cutoff_skin=cutoff_skin)
      qmmm_pot     = set_qmmm_pot(atoms, atoms.params['CrackPos'], mm_pot, qm_pot)

    strain_atoms = fix_edges(atoms)

    print 'Setup dynamics'
#If input_file is crack.xyz the cell has not been thermalized yet.
#Otherwise it will recover temperature from the previous run.
    print 'Attaching trajectories to dynamics'
    trajectory = AtomsWriter(traj_file)
#Only wriates trajectory if the system is in the LOTFDynamicas
#Interpolation 
    atoms.wrap()
    atoms.set_cutoff(3.0)
    atoms.calc_connect()
    print 'Running Crack Simulation'
    RELAXATION = False
    if RELAXATION:
      dynamics = FIRE(atoms)
      dynamics.attach(pass_trajectory_context(trajectory, dynamics), traj_interval, dynamics)
      dynamics.run(fmax=0.1)
    else:
      dynamics = LOTFDynamics(atoms, timestep, extrapolate_steps)
      dynamics.attach(pass_trajectory_context(trajectory, dynamics), traj_interval, dynamics)
      nsteps    = 2000
      dynamics.run(nsteps)
