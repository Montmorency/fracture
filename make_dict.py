import pickle
import ase.units as units

crack_info ={
  'input_file'  : 'frac_cell.xyz',    # crack_slab
  'sim_T'       : 300.0*units.kB, # Simulation temperature
  'nsteps'      : 10000,          # Total number of timesteps to run for
  'timestep'    : 1.0*units.fs,  # Timestep (NB: time base units are not fs!)
  'cutoff_skin' : 2.0*units.Ang, # Amount by which potential cutoff is increased
                                 # for neighbour calculations
  'tip_move_tol'      : 10.0,    # Distance tip has to move before crack
                                 # is taken to be running
  'strain_rate'       : 1e-3*(1.0/units.fs),
  'traj_interval'     : 10,                   # Number of time steps between
  'traj_file'         : 'traj_lotf_2.xyz',    # Trajectory output file in (NetCDF format)
  'restart_traj_file' : 'traj_lotf_2b.xyz',   # Trajectory output file in (NetCDF format)
  'print_interval'    : 10,              # time steps between trajectory prints 10 fs
  'param_file'        : 'Fe_Mendelev.xml', # Filename of XML file containing
                                         # potential parameters
  'mm_init_args'      :'IP EAM_ErcolAd', # Classical potential
  'qm_init_args'      :'TB DFTB',       # Initialisation arguments for QM potential
  'qm_inner_radius'   : 15.0*units.Ang, # Inner hysteretic radius for QM region
  'qm_outer_radius'   : 20.0*units.Ang, # Outer hysteretic radius for QM region
  'extrapolate_steps' : 10,             # Number of steps for predictor-corrector
                                        # interpolation and extrapolation
  'crack_front'       : (1,1,0),
  'crack_direction'   : (-12,12,-2), 
  'cleavage_plane'    : (-1,1,12),
  'symbol'            : 'Fe',
  'width'  : 900.0*units.Ang,      # Width of crack slab
  'height' : 300.00*units.Ang,      # Height of crack slab
  'vacuum' : 20.0*units.Ang,        # Amount of vacuum around slab
  'crack_seed_length'  : 450.0*units.Ang,   # Length of seed crack
  'strain_ramp_length' : 80.0*units.Ang,    # Distance over which strain is ramped up
  'initial_G'  : 4.3*(units.J/units.m**2),  # Initial energy flow to crack tip
  'relax_fmax' : 0.025*units.eV/units.Ang   # Maximum force criteria
}

with open('crack_info.pckl','w') as f:
  pickle.dump(crack_info,f)
