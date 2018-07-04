import pickle
import ase.units as units
import json

crack_info ={
  'input_file'  : 'frac_cell.xyz',    # crack_slab
  'sim_T'       : 0.0*units.kB, # Simulation temperature
  'nsteps'      : 10000,          # Total number of timesteps to run for
  'timestep'    : 0.25*units.fs,  # Timestep (NB: time base units are not fs!)
  'cutoff_skin' : 2.0*units.Ang, # Amount by which potential cutoff is increased
                                 # for neighbour calculations
  'tip_move_tol'      : 10.0,    # Distance tip has to move before crack
                                 # is taken to be running
  'strain_rate'       : 0.0*(1.0/units.fs),
  'traj_interval'     : 5,                   # Number of time steps between
  'traj_file'         : 'crack_traj.xyz',    # Trajectory output file in (NetCDF format)
  'restart_traj_file' : 'crack_traj.xyz',   # Trajectory output file in (NetCDF format)
  'print_interval'    : 1,              # time steps between trajectory prints 10 fs
  'param_file'        : 'PotBH.xml', # Filename of XML file containing
                                         # potential parameters
  'mm_init_args'      : 'IP EAM_ErcolAd', # Classical potential
  'qm_init_args'      : 'VASP DFT',       # Initialisation arguments for QM potential
  'qm_inner_radius'   : 3.5*units.Ang, # Inner hysteretic radius for QM region
  'qm_outer_radius'   : 15.0*units.Ang, # Outer hysteretic radius for QM region
  'extrapolate_steps' : 10,             # Number of steps for predictor-corrector
                                      # interpolation and extrapolation
#Notation (cleavage_plane)[crack_front]
  'cleavage_plane'    : (1,1,1),
  'crack_front'       : (1,1,2),
  'crack_direction'   : (1,-1,0), 
  'symbol'            : 'Fe',
  'width'  : 600.0*units.Ang,      # Width of crack slab
  'height' : 200.00*units.Ang,      # Height of crack slab
  'vacuum' : 20.0*units.Ang,        # Amount of vacuum around slab
  'crack_seed_length'  : 300.0*units.Ang,   # Length of seed crack
  'strain_ramp_length' : 100.0*units.Ang,    # Distance over which strain is ramped up
  'initial_G'  : 4.0*(units.J/units.m**2),  # Initial energy flow to crack tip
  'relax_fmax' : 0.01*units.eV/units.Ang   # Maximum force criteria
}

with open('crack_info.json','w') as f:
  json.dump(crack_info,f,indent=1)
