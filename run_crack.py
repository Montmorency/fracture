import os 
import json
import argparse
import numpy as np
import ase.units as units



from imeall import app

from ase.io import Trajectory
from ase.io.xyz import write_xyz
from ase import Atoms as aseAtoms
from ase.constraints import FixAtoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from quippy           import Atoms
from quippy.potential import Potential
from quippy.io        import AtomsWriter, AtomsReader
from quippy.crack     import get_strain, get_energy_release_rate,\
                             ConstantStrainRate, find_crack_tip_stress_field

from quippy           import set_fortran_indexing

from quippy.potential import ForceMixingPotential
from quippy.lotf      import LOTFDynamics, update_hysteretic_qm_region

from matscipy.socketcalc import VaspClient, SocketCalculator

from ForceMixerCarver import ForceMixingCarvingCalculator
from simulate_crack   import update_qm_region_context, fix_edges, set_qmmm_pot, pass_print_context,\
                             check_if_cracked_context, pass_trajectory_context

from distutils import spawn

set_fortran_indexing(False)
#simulation parameters
extrapolate_steps = 16         # Number of steps for predictor-corrector
                               # interpolation and extrapolation
with open('crack_info.json','r') as f:
  crack_dict = json.load(f)

sim_T       = crack_dict['sim_T'] # Simulation temperature
nsteps      = 6000             # Total number of timesteps to run for
timestep    = 1.0*units.fs    # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang    # Amount by which potential cutoff is increased
                               # for neighbour calculations
tip_move_tol = 12.0            # Distance tip has to move before crack
                               # is taken to be running
strain_rate    = 0.0
traj_interval  = 80             # Number of time steps between interpolations
print_interval = 10             # time steps between trajectory prints 10 fs
param_file     = 'params.xml'   # Filename of XML file containing
                                # potential parameters
mm_init_args = 'IP EAM_ErcolAd do_rescale_r=T r_scale=1.008948312' # Classical potential
qm_init_args = 'TB DFTB'        # Initialisation arguments for QM potential

input_file   = 'crack.xyz'      # crack_slab
traj_file    = 'crack_traj.xyz' # Trajectory output file in (NetCDF format)

try:
  pot_dir      = os.environ['POTDIR']
except:
  sys.exit("PLEASE SET export POTDIR='path/to/potfiles/')")

pot_file     = os.path.join(pot_dir, 'PotBH.xml')

# Restart from last point in trajectory file:
if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input_file",  help="File containing initial configuration. \
                                                   Default is crack.xyz.", default='crack.xyz')
  parser.add_argument("-o", "--output_file", help="File trajectory is written to.", default='./crack_traj.xyz')
  parser.add_argument("-p", "--pot_file",    help="Potential is set.", default=os.path.join(pot_dir,'PotBH.xml'))
  parser.add_argument("-r", "--restart",     help="If false thermalizes atoms and fixes boundary conditions,\
                                                   frame in trajectory.", action="store_true")
  parser.add_argument("-l", "--lotf", help="If true do a LOTF simulation.", action="store_true")
  parser.add_argument("-s", "--socket", help="If true do a verlet simulation with the ForceMixer.", action="store_true")

  args = parser.parse_args()

  if args.restart:
    print 'Restarting job'

  pot_file   = args.pot_file
  if args.output_file:
    traj_file = args.output_file

#Need to add an arg parser here
  if args.socket:
    print 'Initialize Potentials'
    print 'Read last MD snapshot'
    print 'Reading from ', args.input_file

    mpirun = spawn.find_executable('mpirun')
    vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'
#mm_pot
    POT_DIR = os.path.join(app.root_path, 'potentials')
    eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
    r_scale = 1.00894848312
    mm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)

    atoms = AtomsReader(args.input_file)[-1]
    strain_atoms = fix_edges(atoms)
    crack_pos = atoms.info['CrackPos']
    x, y, z = atoms.positions.T
    radius1 = np.sqrt((x - crack_pos[0])**2 + (y-crack_pos[1])**2 + (z-crack_pos[2])**2)

    qm_radius = 3.0
    buff = 9.0

    qm_region_mask = (radius1 < qm_radius)
    qm_buffer_mask = (radius1 < qm_radius + buff)

    print ("\nNumber of atoms in qm region of %.1f" % qm_radius +
                                    "A : %i" % np.count_nonzero(qm_region_mask))

    print ("together with the buffer of %.1f" % (qm_radius + buff ) +
                                    "A %i" % np.count_nonzero(qm_buffer_mask))

    vasp_args = dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
                     kpts=[1, 1, 6], kpar=1, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-4,
                     nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0, encut=400,
                     maxmix=30, voskown=0, ismear=1, sigma=0.1, isym=0)
    procs = 96
# need to have procs % n_par == 0
    n_par = 1
    if procs <= 8:
        n_par = procs
    else:
        for _npar in range(2, int(np.sqrt(1.*procs))):
            if procs % int(_npar) == 0:
                n_par = procs // int(_npar)
    vasp_client = VaspClient(client_id=0, npj=procs, ppn=1,
                             exe=vasp, mpirun=mpirun, parmode='mpi',
                             ibrion=13, nsw=1000000,
                             npar=n_par, **vasp_args)
    qm_pot = SocketCalculator(vasp_client)
    qmmm_pot = ForceMixingCarvingCalculator(atoms, qm_region_mask,
                                            mm_pot, qm_pot,
                                            buffer_width=buff,
                                            pbc_type=[False, False, True])
    atoms.set_calculator(qmmm_pot)
#Otherwise it will recover temperature from the previous run.
#Use same random seed so that initialisations are deterministic.
    if not args.restart: 
        print 'Thermalizing atoms'
        np.random.seed(42)
        MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)
    dynamics = VelocityVerlet(atoms, timestep)

    def print_context(ats=atoms, dyn=dynamics):
        print 'steps, T', dyn.nsteps, ats.get_kinetic_energy()/(1.5*units.kB*len(ats))
        print 'G', get_energy_release_rate(ats)/(units.J/units.m**2)
        print 'strain', get_strain(ats)
    dynamics.attach(print_context, interval=8)
    print 'Running Crack Simulation'
    traj = Trajectory('run_crack.traj', 'w', atoms)
    dynamics.attach(traj.write, interval=32)
    dynamics.run(nsteps)
    print 'Crack Simulation Finished'

  else:
    mm_pot = Potential(mm_init_args, param_filename=pot_file, cutoff_skin=cutoff_skin)
    atoms = AtomsReader(args.input_file)[-1]
    atoms.set_calculator(mm_pot)
    strain_atoms = fix_edges(atoms)
    current_crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot)
    print 'Current Crack Position: ', current_crack_pos
    if not args.restart: 
      print 'Thermalizing Atoms'
      MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)
    dynamics = VelocityVerlet(atoms, timestep)

    def print_context(ats=atoms, dyn=dynamics):
        print 'steps, T', dyn.nsteps, ats.get_kinetic_energy()/(1.5*units.kB*len(ats))
        print 'G', get_energy_release_rate(ats)/(units.J/units.m**2)
        print 'strain', get_strain(ats)
    dynamics.attach(print_context, interval=8)
    print 'Running Crack Simulation'
    #traj = Trajectory('run_crack_IV.traj', 'w', atoms)
    #dynamics.attach(traj.write, interval=32)
    def write_slab(a=atoms):
        write_xyz('crack_traj.xyz', a, append=True)
    dyanmics.attach(traj.write, interval=8)
    dynamics.run(nsteps)
    print 'Crack Simulation Finished'
