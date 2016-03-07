from cStringIO import StringIO
import sys
import pickle

import ase.units as units
from   ase.lattice       import bulk
from   ase.lattice.cubic import Diamond, BodyCenteredCubic
from   ase.constraints   import FixAtoms

from   quippy import set_fortran_indexing
from   quippy.potential  import Potential, Minim
from   quippy.elasticity import youngs_modulus, poisson_ratio, rayleigh_wave_speed, AtomResolvedStressField
from   quippy.io    import write
from   quippy.crack import (print_crack_system,
                            G_to_strain,
                            thin_strip_displacement_y,
                            find_crack_tip_stress_field)

from quippy import bcc, calc_nye_tensor, Atoms, supercell
from qlab import view, aux_property_coloring

import numpy as np

set_fortran_indexing(False)
#
# Crack Object for rapid generation of prototypical crack simulation
# in an ideal crystal. This should allow for quick generation of
# different crack geometries, strain rates, etc. The crack_params
# dictionary should take care of all the parameters.
#
crack_params = {
	'input_file'  : 'crack.xyz',    # crack_slab
	'sim_T'       : 300.0*units.kB, # Simulation temperature
	'nsteps'      : 10000,          # Total number of timesteps to run for
	'timestep'    : 1.0*units.fs,  # Timestep (NB: time base units are not fs!)
	'cutoff_skin' : 2.0*units.Ang, # Amount by which potential cutoff is increased
                                 # for neighbour calculations
	'tip_move_tol'      : 10.0,    # Distance tip has to move before crack
                                 # is taken to be running
	'strain_rate'       : 1e-5*(1.0/units.fs),
	'traj_interval'     : 10,                   # Number of time steps between
	'traj_file'         : 'traj_lotf_2.xyz',    # Trajectory output file in (NetCDF format)
	'restart_traj_file' : 'traj_lotf_2b.xyz',   # Trajectory output file in (NetCDF format)
	'print_interval'    : 10,              # time steps between trajectory prints 10 fs
	'param_file'        : 'iron_mish.xml', # Filename of XML file containing
                                         # potential parameters
	'mm_init_args'      :'IP EAM_ErcolAd', # Classical potential
	'qm_init_args'      :'TB DFTB',      # Initialisation arguments for QM potential
	'qm_inner_radius'   : 15.0*units.Ang, # Inner hysteretic radius for QM region
	'qm_outer_radius'   : 20.0*units.Ang, # Outer hysteretic radius for QM region
	'extrapolate_steps' : 10,      	 	    # Number of steps for predictor-corrector
  				                              # interpolation and extrapolation
#	'cleavage_plane'    : (1,0,0),
#	'crack_front'       : (0,1,0),
#	'crack_direction'   : (0,0,1),
	'cleavage_plane'    : (1,-1,0),
	'crack_front'       : (1,1,0),
	'crack_direction'   : (0,0,1),

	'symbol'            : 'Fe',

	'width'  : 300.0*units.Ang,        # Width of crack slab
	'height' : 150.0*units.Ang,        # Height of crack slab
	'vacuum' : 100.0*units.Ang,        # Amount of vacuum around slab

	'crack_seed_length'  : 60.0*units.Ang,    # Length of seed crack
	'strain_ramp_length' : 60.0*units.Ang,    # Distance over which strain is ramped up
	#'initial_G'  : 5.0*(units.J/units.m**2), # Initial energy flow to crack tip
	'initial_G'  : 11.4*(units.J/units.m**2), # Initial energy flow to crack tip
	'relax_fmax' : 0.025*units.eV/units.Ang   # Maximum force criteria
}

#Return dictionary of default params
#these can be edited and then pickled.
#with gen_inputfile
def init_dict():
	return crack_params

def gen_inputfile(crack_dict, filename):
	f = open(filename.join('.pckl'), 'w')
	g = open(filename.join('txt', 'w'))
	for key, value in crack_dict.items():
		print >> g, key, value
	pickle.dump(f, crack_dict)
	f.close()
	g.close()

class Capturing(list):
	def __enter__(self):
		self._stdout = sys.stdout
		sys.stdout = self._stringio = StringIO()
		return self
	def __exit__(self, *args):
		self.extend(self._stringio.getvalue().splitlines())
		sys.stdout = self._stdout

class CrackCell(object):
	def __init__(self, **kwargs):
		self.symbol = 'Fe'
		self.width  = 300.0*units.Ang        # Width of crack slab
		self.height = 150.0*units.Ang        # Height of crack slab
		self.vacuum = 100.0*units.Ang        # Amount of vacuum around slab
		self.crack_seed_length  = 60.0*units.Ang    # Length of seed crack
		self.strain_ramp_length = 45.0*units.Ang    # Distance over which strain is ramped up
		self.initial_G  = 5.0*(units.J/units.m**2)  # Initial energy flow to crack tip
		self.relax_fmax = 0.025*units.eV/units.Ang  # Maximum force criteria
		self.param_file = 'params.xml'              # XML file containing
                                   # interatomic potential parameters
		self.mm_init_args    = 'IP SW' # Initialisation arguments
		                               # for the classical potential
		self.output_file     = 'crack.xyz'
		self.cleavage_plane  = (-2, 1, 1)
		self.crack_direction = (0, 1, -1)
		self.crack_front     = (1,1,1)
		self.a0              = 2.83
		self.nx              = 1
		self.ny              = 1
		self.c               = []

#Initialize the crack object with a dictionary of the relevant parameters
		for key in kwargs:
			setattr(self, key, kwargs[key])

	def __repr__(self):
		with Capturing() as output:
			print_crack_system(self.crack_direction, self.cleavage_plane, self.crack_front)
		return '\n'.join(output)
	
	def calculate_c(self,size=(1,1,1)):
		if self.symbol=='Si':
			bulk = Diamond(directions = [self.crack_direction, self.cleavage_plane, self.crack_front],
		  							 size = size, symbol=self.symbol, pbc=True, latticeconstant= self.a0)

		elif self.symbol=='Fe':
			bulk = BodyCenteredCubic(directions=[self.crack_direction, self.cleavage_plane, self.crack_front],
                               size=size, symbol='Fe', pbc=(1,1,1),
                               latticeconstant=2.83)

		pot     = Potential(self.mm_init_args, param_filename=self.param_file)
		bulk.set_calculator(pot)
		self.c  = pot.get_elastic_constants(bulk)
		return

	def build_unit_slab(self, size=(1,1,1)):
		if self.symbol=='Si':
			unit_slab = Diamond(directions = [self.crack_direction, self.cleavage_plane, self.crack_front],
													size = size, symbol=self.symbol, pbc=True, latticeconstant= self.a0)
		elif self.symbol=='Fe':
			unit_slab = BodyCenteredCubic(directions=[self.crack_direction, self.cleavage_plane, self.crack_front],
                               size=size, symbol='Fe', pbc=(1,1,1),
                               latticeconstant=2.83)
#does this work for more than 2 atoms?
		unit_slab.positions[:, 1] += (unit_slab.positions[1, 1]-unit_slab.positions[0, 1])/2.0
		unit_slab.set_scaled_positions(unit_slab.get_scaled_positions())
		print self.mm_init_args, self.param_file
		pot     = Potential(self.mm_init_args, param_filename=self.param_file)
		unit_slab.set_calculator(pot)
		return unit_slab

	def build_surface(self, size=(1,1,1)):
		if self.symbol=='Si':
			unit_slab = Diamond(directions = [self.crack_direction, self.cleavage_plane, self.crack_front],
													size = size, symbol=self.symbol, pbc=True, latticeconstant= self.a0)
		elif self.symbol=='Fe':
			unit_slab = BodyCenteredCubic(directions=[self.crack_direction, self.cleavage_plane, self.crack_front],
                               size=size, symbol='Fe', pbc=(1,1,1),
                               latticeconstant=self.a0)
#does this work for more than 2 atoms?
		unit_slab.positions[:, 1] += (unit_slab.positions[1, 1]-unit_slab.positions[0, 1])/2.0
		unit_slab.set_scaled_positions(unit_slab.get_scaled_positions())
		surface = unit_slab.copy()
		surface.center(vacuum=self.vacuum, axis=1)
		pot     = Potential(self.mm_init_args, param_filename=self.param_file)
		surface.set_calculator(pot)
		return surface

	def build_crack_cell(self, unit_slab):
#return a copy of the surface and the crack_cell
#with the constraints on the Atoms.
		self.nx = int(self.width/unit_slab.get_cell()[0,0])
		self.ny = int(self.height/unit_slab.get_cell()[1,1])
		if self.ny %2 == 1:
			self.ny +=1
		print 'number of unit cells', self.nx, self.ny
		crack_slab = unit_slab*(self.nx, self.ny,1)
		crack_slab.center(self.vacuum, axis=0)
		crack_slab.center(self.vacuum, axis=1)
		write('crack_slab_orig.xyz',crack_slab)
		crack_slab.positions[:, 0] -= crack_slab.positions[:, 0].mean()
		crack_slab.positions[:, 1] -= crack_slab.positions[:, 1].mean()
		orig_width  = (crack_slab.positions[:, 0].max() -
		               crack_slab.positions[:, 0].min())
		orig_height = (crack_slab.positions[:, 1].max() -
		               crack_slab.positions[:, 1].min())

		print(('Made slab with %d atoms, original width and height: %.1f x %.1f A^2' %
       (len(crack_slab), orig_width, orig_height)))

		top    = crack_slab.positions[:,1].max()
		bottom = crack_slab.positions[:,1].min()
		left   = crack_slab.positions[:,0].min()
		right  = crack_slab.positions[:,0].max()
#Fix top and bottom layer of atoms:
		fixed_mask = ((abs(crack_slab.positions[:, 1] - top) < 1.0) |
                  (abs(crack_slab.positions[:, 1] - bottom) < 1.0))
		const = FixAtoms(mask=fixed_mask)
		crack_slab.set_constraint(const)
#print('Fixed %d atoms\n' % fixed_mask.sum())
#Strain is a dimensionless ratio require elastic tensor
#to translate the initial energy flow to the tip into a strain
		E  = youngs_modulus(self.c, self.cleavage_plane)
		nu = poisson_ratio(self.c, self.cleavage_plane, self.crack_direction)
		strain = G_to_strain(self.initial_G, E, nu, orig_height)
		seed = left + self.crack_seed_length
		tip  = left + self.crack_seed_length + self.strain_ramp_length
		xpos = crack_slab.positions[:,0]
		ypos = crack_slab.positions[:,1]
		crack_slab.positions[:,1] += thin_strip_displacement_y(xpos, ypos, strain, seed, tip)
		print('Applied initial load: strain=%.4f, G=%.2f J/m^2' %
	       (strain, self.initial_G / (units.J / units.m**2)))
		pot      = Potential(self.mm_init_args, param_filename=self.param_file)
		crack_slab.set_calculator(pot)
		slab_opt = Minim(crack_slab, relax_positions=True, relax_cell=False)
		slab_opt.run(fmax=self.relax_fmax)
		return crack_slab

	def calc_nye_tensor(self, x):
		x.set_cutoff(3)
		x.calc_connect()
		b = bcc(2.83, 26)
		x0 = supercell(b, self.nx, self.ny, 1)
		x0.set_cutoff(3)
		x0.calc_connect()
		alpha = calc_nye_tensor(x,x0,3,3,x.n)
		a = np.zeros([3,3])
		for i in range(3):
			for j in range(3):
				a[i,j] = sum(alpha[i,j,:])

		return alpha, a

if __name__ == '__main__':
#	print crack_params
	crack      = CrackCell(**crack_params)
	unit_slab  = crack.build_unit_slab()
	crack.calculate_c()
	surface    = crack.build_surface()
	print surface.get_cell()
	E_surf = surface.get_potential_energy()
	E_bulk = unit_slab.get_potential_energy()/len(unit_slab)
	area = (surface.get_cell()[0,0]*surface.get_cell()[2,2])
	print E_surf, E_bulk
	gamma = (E_surf - E_bulk*len(surface))/(2.0*area)
	crack_slab = crack.build_crack_cell(unit_slab)
	print('Surface energy of %s surface %.4f J/m^2\n' %
      (crack.cleavage_plane, gamma / (units.J / units.m ** 2)))
	print crack_slab
	write('crack_slab.xyz', crack_slab)

#	alpha, burg = crack.calc_nye_tensor(Atoms)
#	print burg.round(4)
#	v = view(Atoms('crack_slab.xyz'))
#	aux_property_coloring(alpha[1,2,:])
#	v.show()

