"""
Module containing CrackCell object.
"""
from cStringIO import StringIO
import sys
import json

import ase.units as units
from ase.lattice import bulk
from ase.constraints import FixAtoms
from ase.lattice.cubic import Diamond, BodyCenteredCubic
from ase.optimize import FIRE
#from  ase.optimize      import LBFGS
#preconditioned LBFGS is !much! faster for the crack slab
#relaxations tested so far.
from ase.optimize.precon  import PreconLBFGS, PreconFIRE
from quippy import set_fortran_indexing
from quippy.elasticity import youngs_modulus, poisson_ratio, rayleigh_wave_speed, AtomResolvedStressField
from quippy.io import write
from quippy.crack import print_crack_system, G_to_strain,\
                            thin_strip_displacement_y,\
                            find_crack_tip_stress_field

from quippy import Potential
from quippy import bcc, calc_nye_tensor, Atoms, supercell

import os
import numpy as np

set_fortran_indexing(False)
default_crack_params = {
  'input_file'  : 'crack.xyz',    # crack_slab
  'sim_T'       : 300.0*units.kB, # Simulation temperature
  'nsteps'      : 10000,          # Total number of timesteps to run for
  'timestep'    : 1.0*units.fs,   # Timestep (NB: time base units are not fs!)
  'cutoff_skin' : 2.0*units.Ang,  # Amount by which potential cutoff is increased
                                  # for neighbour calculations.
  'tip_move_tol'      : 10.0,     # Distance tip has to move before crack
                                  # is taken to be running.
  'strain_rate'       : 1e-7*(1.0/units.fs),
  'traj_interval'     : 10,                 # Number of time steps between interpolations
  'traj_file'         : 'traj_lotf_2.xyz',  # Trajectory output file in (NetCDF format)
  'restart_traj_file' : 'traj_lotf_2b.xyz', # Trajectory output file in (NetCDF format)
  'print_interval'    :  20,              # time steps between trajectory prints n fs
  'param_file'        : 'params.xml', # Filename of XML file containing
                                         # potential parameters
  #'mm_init_args'      :'IP EAM_ErcolAd', # Classical potential
  'mm_init_args'      :'IP SW', # Classical potential
  'qm_init_args'      :'TB DFTB',       # Initialisation arguments for QM potential
  'qm_inner_radius'   : 18.0*units.Ang, # Inner hysteretic radius for QM region
  'qm_outer_radius'   : 21.0*units.Ang, # Outer hysteretic radius for QM region
  'extrapolate_steps' : 10,               # Number of steps for predictor-corrector
                                        # interpolation and extrapolation
  'cleavage_plane'    : (1,-1,0),
  'crack_front'       : (1,1,0),
  'crack_direction'   : (0,0,1),
  'symbol'            : 'Si',
  'width'  : 900.0*units.Ang,             # Width of crack slab
  'height' : 300.0*units.Ang,             # Height of crack slab
  'vacuum' : 25.0*units.Ang,              # Amount of vacuum around slab
  'crack_seed_length'  : 250.0*units.Ang, # Length of seed crack
  'strain_ramp_length' : 50.0*units.Ang,  # Distance over which strain is ramped up
  'initial_G'  : 4.5*(units.J/units.m**2), # Initial energy flow to crack tip
  'relax_fmax' : 0.025*units.eV/units.Ang   # Maximum force criteria
}

# Return dictionary of default params
# these can be edited and then pickled.
# with gen_inputfile
def init_dict(**kwargs):
  for key in kwargs:
    default_crack_params[key] = kwargs[key]
  return default_crack_params

#Capturing
class Capturing(list):
  def __enter__(self):
    self._stdout = sys.stdout
    sys.stdout = self._stringio = StringIO()
    return self
  def __exit__(self, *args):
    self.extend(self._stringio.getvalue().splitlines())
    sys.stdout = self._stdout

class CrackCell(object):
  """
  CrackCell defines the parameters
  for an ideal mode 1 crack geometry. This is intended 
  to allow for the quick generation of different crack
  geometries, strain rates, etc. It also encapsulates
  the methods to initialize the seed crack and relax the 
  initial configuration.
  """
  def __init__(self, cij=[], **kwargs):
    self.symbol = 'Fe'
    self.width  = 500.0*units.Ang        # Width of crack slab
    self.height = 166.0*units.Ang        # Height of crack slab
    self.vacuum = 100.0*units.Ang        # Amount of vacuum around slab
    self.crack_seed_length  = 250.0*units.Ang    # Length of seed crack
    self.strain_ramp_length = 60.0*units.Ang    # Distance over which strain is ramped up
    self.initial_G  = 4.5*(units.J/units.m**2)  # Initial energy flow to crack tip
    self.relax_fmax = 0.020*units.eV/units.Ang  # Maximum force criteria
    self.param_file = 'params.xml'              # XML file containing
                                   # interatomic potential parameters
    self.mm_init_args    = 'IP SW' # Initialisation arguments
                                   # for the classical potential
    self.output_file     = 'crack.xyz'
    self.cleavage_plane  = (1, 1, 1)
    self.crack_direction = (1, 1, -2)
    self.crack_front     = (1,-1,0)
    self.a0              = 2.8293
    self.nx              = 1
    self.ny              = 1
    self.nz              = 1
#These are derived properties cell must be initialized before
#they are calculated
    self.cij             = cij
    self.E               = 0.0
    self.nu              = 0.0
    self.orig_height     = 0.0
    self.orig_width      = 0.0
    self.strain          = 0.0

#Initialize the crack object with a dictionary of the relevant parameters
    for key in kwargs:
      setattr(self, key, kwargs[key])

  def __repr__(self):
    """
    Print the CrackCell geometry in canonical form.
    """
    with Capturing() as output:
      print_crack_system(self.crack_direction, self.cleavage_plane, self.crack_front)
    return '\n'.join(output)
  
  def calculate_c(self, pot, size=(1,1,1)):
    """
    Calculate the elastic constants for the underlying crack unit cell.
    """
    if self.symbol=='Si':
      at_bulk = bulk('Si', a=self.a0, cubic=True)
    elif self.symbol=='Fe':
      at_bulk = BodyCenteredCubic(directions=[self.crack_direction, self.cleavage_plane, self.crack_front],
                               size=size, symbol='Fe', pbc=(1,1,1),
                               latticeconstant=self.a0)
    at_bulk.set_calculator(pot)
    self.cij  = pot.get_elastic_constants(at_bulk)
    print((self.cij / units.GPa).round(2))
    return

  def build_unit_slab(self, size=(1,1,1)):
    """
    Initialize the unit slab for the given crack geometry.
    Currently supports silicon and Fe.
    """
    if self.symbol=='Si':
      unit_slab = Diamond(directions = [self.crack_direction, self.cleavage_plane, self.crack_front],
                          size = size, symbol=self.symbol, pbc=True, latticeconstant= self.a0)
    elif self.symbol=='Fe':
      unit_slab = BodyCenteredCubic(directions=[self.crack_direction, self.cleavage_plane, self.crack_front],
                               size=size, symbol='Fe', pbc=(1,1,1),
                               latticeconstant=self.a0)
    print 'Number atoms in unit slab', len(unit_slab)
    unit_slab.info['adsorbate_info'] = None
    unit_slab.positions[:, 1] += (unit_slab.positions[1, 1]-unit_slab.positions[0, 1])/2.0
    unit_slab.set_scaled_positions(unit_slab.get_scaled_positions())
    pot_dir = os.environ['POTDIR']
    pot_file = os.path.join(pot_dir, self.param_file)
    print pot_file
    print self.mm_init_args
    #pot = Potential(self.mm_init_args, param_filename = pot_file)
    #pot = Potential("IP EAM_ErcolAd do_rescale_r=T r_scale=1.00894848312" , param_filename = "/Users/lambert/pymodules/imeall/imeall/potentials/PotBH.xml")
    #unit_slab.set_calculator(pot)
    return unit_slab

  def build_surface(self,pot, size=(1,1,1)):
    """
    Create an equivalent surface unit cell for the fracture system.
    """
    if self.symbol=='Si':
      unit_slab = Diamond(directions = [self.crack_direction, self.cleavage_plane, self.crack_front],
                          size = size, symbol=self.symbol, pbc=True, latticeconstant= self.a0)
    elif self.symbol=='Fe':
      unit_slab = BodyCenteredCubic(directions=[self.crack_direction, self.cleavage_plane, self.crack_front],
                               size=size, symbol='Fe', pbc=(1,1,1),
                               latticeconstant=self.a0)
# Does this work for more than 2 atoms?
    unit_slab.info['adsorbate_info'] = None
    unit_slab.positions[:, 1] += (unit_slab.positions[1, 1]-unit_slab.positions[0, 1])/2.0
    unit_slab.set_scaled_positions(unit_slab.get_scaled_positions())
    surface = unit_slab.copy()
    surface.center(vacuum=self.vacuum, axis=1)
    surface.set_calculator(pot)
    return surface

  def build_crack_cell(self, unit_slab, pot, fix_dist=1.0):
    """
    Creates a CrackCell of the desired size and orientation
    with the top and bottom layer of atoms constrained.
    Returns: :class:`Atoms` - the relaxed crack_cell.
    """
    self.nx = int(self.width/unit_slab.get_cell()[0,0])
    self.ny = int(self.height/unit_slab.get_cell()[1,1])
    if self.ny %2 == 1:
      self.ny +=1
    #self.ny = 1
    print 'number of unit cells', self.nx, self.ny
    crack_slab = unit_slab*(self.nx, self.ny,self.nz)
    write('ref_slab.xyz', crack_slab)
    crack_slab.center(self.vacuum, axis=0)
    crack_slab.center(self.vacuum, axis=1)
    crack_slab.positions[:, 0] -= crack_slab.positions[:, 0].mean()
    crack_slab.positions[:, 1] -= crack_slab.positions[:, 1].mean()
    self.orig_width  = (crack_slab.positions[:, 0].max() -
                   crack_slab.positions[:, 0].min())
    self.orig_height = (crack_slab.positions[:, 1].max() -
                   crack_slab.positions[:, 1].min())
    print(('Made slab with %d atoms, original width and height: %.1f x %.1f A^2' %
       (len(crack_slab), self.orig_width, self.orig_height)))
    top    = crack_slab.positions[:,1].max()
    bottom = crack_slab.positions[:,1].min()
    left   = crack_slab.positions[:,0].min()
    right  = crack_slab.positions[:,0].max()
#Constrain top and bottom layer of atoms:
    fixed_mask = ((abs(crack_slab.positions[:, 1] - top) < 1.0) |
                  (abs(crack_slab.positions[:, 1] - bottom) < 1.0))
    const = FixAtoms(mask=fixed_mask)
    crack_slab.set_constraint(const)

#Strain is a dimensionless ratio we require elastic tensor
#to translate the initial energy flow to the tip into a strain.

    self.E  = youngs_modulus(self.cij, self.cleavage_plane)
    self.nu = poisson_ratio(self.cij, self.cleavage_plane, self.crack_direction)
    self.strain = G_to_strain(self.initial_G, self.E, self.nu, self.orig_height)
    seed = left + self.crack_seed_length
    tip  = left + self.crack_seed_length + self.strain_ramp_length
    xpos = crack_slab.positions[:,0]
    ypos = crack_slab.positions[:,1]
    crack_slab.positions[:,1] += thin_strip_displacement_y(xpos, ypos, self.strain, seed, tip)
    write('crack_slab_init.xyz', crack_slab)

    print('Applied initial load: strain=%.4f, G=%.2f J/m^2' %
         (self.strain, self.initial_G / (units.J / units.m**2)))

    pot_dir  = os.environ['POTDIR']
    crack_slab.set_calculator(pot)
#Initial crack pos
    print 'crack_pos before relaxations'
    print  find_crack_tip_stress_field(crack_slab, calc=mm_pot)
    print('Relaxing slab...')
    slab_opt = PreconFIRE(crack_slab) 
    slab_opt.run(fmax=self.relax_fmax)
    return crack_slab

  def calc_nye_tensor(self, x):
    x.set_cutoff(3)
    x.calc_connect()
    b = bcc(2.82893, 26)
    x0 = supercell(b, self.nx, self.ny, 1)
    x0.set_cutoff(3)
    x0.calc_connect()
    alpha = calc_nye_tensor(x,x0,3,3,x.n)
    a = np.zeros([3,3])
    for i in range(3):
      for j in range(3):
        a[i,j] = sum(alpha[i,j,:])
    return alpha, a

  def write_crack_cell(self, crack_slab, mm_pot):
    print crack_slab.get_cell()
    crack_slab.info['nneightol']       = 1.30 # set nearest neighbour tolerance
    crack_slab.info['LatticeConstant'] = self.a0
    crack_slab.info['C11'] = self.cij[0, 0]
    crack_slab.info['C12'] = self.cij[0, 1]
    crack_slab.info['C44'] = self.cij[3, 3]
    crack_slab.info['YoungsModulus']   = self.E
    crack_slab.info['PoissonRatio_yx'] = self.nu
    crack_slab.info['G']               = self.initial_G
    crack_slab.info['OrigWidth']       = self.orig_width
    crack_slab.info['OrigHeight']      = self.orig_height
    crack_slab.info['CrackDirection']  = self.crack_direction
    crack_slab.info['CleavagePlane']   = self.cleavage_plane
    crack_slab.info['CrackFront']      = self.crack_front
    crack_slab.info['strain']          = self.strain
    #Initial guess for crack_pos
    crack_slab.info['CrackPos']        = np.array([1.5, 1.5, 1.5])
    crack_pos = find_crack_tip_stress_field(crack_slab, calc=mm_pot)
    print 'Found crack tip at position %s' % crack_pos
    crack_slab.info['CrackPos']        = crack_pos
    crack_slab.info['is_cracked']      = False
    write('crack.xyz', crack_slab)

if __name__ == '__main__':
  #load CrackCell Dictionary from 'crack_info.json' file must be present
  try:
    with open('crack_info.json', 'r') as f:
        crack_info = json.load(f)
    print 'Initializing crack_cell from Info File.'
    crack = CrackCell(**crack_info)
    f.close()
  except IOError:
    print 'no crack dictionary found.'
    sys.exit()

  Grain_Boundary = False

  if not Grain_Boundary:
    unit_slab  = crack.build_unit_slab()

  pot_dir  = os.environ['POTDIR']
  pot_file = os.path.join(pot_dir, crack_info['param_file'])
  mm_pot   = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale=1.00894848312', param_filename=pot_file, cutoff_skin=2.0)

# gb_frac.py will generate the frac_cell.xyz file which contains 
# a crack_cell:
  if Grain_Boundary:
    unit_slab   = Atoms('frac_cell.xyz')
    unit_slab.set_calculator(mm_pot)

#calculate the elasticity tensor:
  crack.calculate_c(mm_pot)
  surface    = crack.build_surface(mm_pot)
  E_surf = surface.get_potential_energy()
  bulk = bcc(2.82893)
  bulk.set_atoms(26)
  bulk.info['adsorbate_info']=None
  bulk.set_calculator(mm_pot)
  E_bulk = bulk.get_potential_energy()/len(bulk)
  area = (surface.get_cell()[0,0]*surface.get_cell()[2,2])
  print E_surf, E_bulk
  gamma = (E_surf - E_bulk*len(surface))/(2.0*area)
  print('Surface energy of %s surface %.4f J/m^2\n' %
       (crack.cleavage_plane, gamma/(units.J/units.m**2)))
  crack_slab = crack.build_crack_cell(unit_slab, mm_pot)
  crack.write_crack_cell(crack_slab, mm_pot)
