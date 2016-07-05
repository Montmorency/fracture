import os
import json
import argparse

import ase.units as units
from   ase.constraints             import FixAtoms
from   ase.lattice.cubic import Diamond, BodyCenteredCubic
from   ase.md.verlet               import VelocityVerlet
from   ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from   quippy            import Atoms, Potential, calc_nye_tensor, set_fortran_indexing
from   quippy.io         import AtomsReader, AtomsWriter
from   quippy.crack      import thin_strip_displacement_y, G_to_strain, get_strain, strain_to_G
from   quippy.elasticity import youngs_modulus, poisson_ratio, rayleigh_wave_speed, AtomResolvedStressField
from   quippy.lotf       import LOTFDynamics, update_hysteretic_qm_region
from   quippy.potential  import ForceMixingPotential, Potential

from   simulate_crack    import update_qm_region_context, set_qmmm_pot, pass_print_context,\
			   							  	   	  check_if_cracked_context, pass_trajectory_context
import numpy as np

set_fortran_indexing = False

# Function to identify the core region: 
def fix_edges_defect(atoms, strain_at=False):
  orig_height    = atoms.info['OrigHeight']
  top    = atoms.positions[:,1].max()
  bottom = atoms.positions[:,1].min()
  left   = atoms.positions[:,0].min()
  right  = atoms.positions[:,0].max()
  fixed_mask = ((abs(atoms.positions[:, 1] - top) < 2.0) |
               (abs(atoms.positions[:, 1] - bottom) < 2.0))
  fix_atoms  = FixAtoms(mask=fixed_mask)
  print('Fixed %d atoms\n' % fixed_mask.sum()) # Print the number of fixed atoms
  if strain_at:
    strain_atoms = ConstantStrainRate(orig_height, strain_rate*timestep)
    atoms.set_constraint([fix_atoms, strain_atoms])
    return strain_atoms
  else:
    atoms.set_constraint([fix_atoms])

def defect_json(**json_args):
  try:
    with open('defect.json','r') as f:
      j_object = json.load(f)
  except IOError:
    j_object = {}
  for key, value in json_args.items():
    j_object[key] = value
  with open('defect.json','w') as f:
    json.dump(j_object, f, indent=2)

def pp_nye_tensor(cl):
#post processing routine to append nye tensor information
#to atoms object.
#Load reference slab and calculate connectivity
    ref_slab = Atoms('./ref_slab.xyz')
    cl.set_cutoff(3.0)
    cl.calc_connect()
    ref_slab.set_cutoff(3.0)
    ref_slab.calc_connect()
    core     = np.zeros(3)
#To save time it is possible to only
#calculate the nye tensor for a 'core' region
#determined by cutting out a cluster:
#    core[:] = atoms.params['core']
#    cl = atoms.select( (atoms.pos[1,:]-core[1])**2 + (atoms.pos[2,:]-core[2])**2 < rr**2)
#Append screw and edge information
#Probably should just pass an array at this stage?
    alpha     = calc_nye_tensor(cl,ref_slab,3,3,cl.n)
    cl.screw  = alpha[2,2,:]
    cl.edgex  = alpha[2,0,:]
    cl.edgey  = alpha[2,1,:]
#Update quantum region according to the
#position of the dislocation core:
    sum = 0
    c=np.array([0.,0.,0.])
    for i in range(cl.n):
        sum = sum + cl.screw[i]
        c[0] = c[0] + cl.pos[0,i]*cl.screw[i]
        c[1] = c[1] + cl.pos[1,i]*cl.screw[i]
    c[0] = c[0]/sum
    c[1] = c[1]/sum
    c[2] = cl.lattice[2,2]/2.
    core[:] = c.copy()
    return


# Initial Parameters:
sim_T      = 300.0*units.kB
input_file = 's111111.xyz'
traj_file  = 's111111_traj.xyz'
timestep       =  5.0*units.fs #Timestep (NB: time base units are not fs!)
print_interval = 20
initial_G      = 0.0*(units.J/units.m**2) #Initial energy flow to crack tip or strain energy of slab
nsteps         = 10000  # Total number of timesteps to run for

#Create unit cell with the orientation:

x = [1,1,-2]
y = [-1,1,0]
z = [1,1,1]
screw_slab_unit = BodyCenteredCubic(directions = [x, y, z], size=(1,1,1),
                                    symbol='Fe', pbc=(1,1,1), latticeconstant=2.83)

POT_DIR  = '/Users/lambert/pymodules/imeall/imeall/potentials'
eam_pot = os.path.join(POT_DIR, 'Fe_Mendelev.xml')
r_scale = 1.00894848312
pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)


if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-rd", "--run_dyn", action='store_true')
  parser.add_argument("-cn", "--calc_nye", action='store_true')
  parser.add_argument("-inp", "--input_file", required=True)
  args = parser.parse_args()
  run_dyn     = args.run_dyn
  calc_nye    = args.calc_nye
  input_file  = args.input_file
  print args.run_dyn, args.calc_nye, input_file

  cij = pot.get_elastic_constants(screw_slab_unit)
  print ((cij / units.GPa).round(2))
  E   = youngs_modulus(cij, y)
  nu  = poisson_ratio(cij, y, x)
  print 'Youngs Modulus: ', E/units.GPa,'Poisson Ratio: ', nu
  print 'Effective elastic modulus E: ', E/(1.-nu**2)
  
  #Load atoms and set potential
  defect = Atoms('{0}.xyz'.format(input_file))
  defect.set_calculator(pot)
  top    = defect.positions[:,1].max()
  bottom = defect.positions[:,1].min()
  left   = defect.positions[:,0].min()
  right  = defect.positions[:,0].max()
  
  orig_height = (defect.positions[:,1].max()-defect.positions[:,1].min())
  defect.info['YoungsModulus']   = E
  defect.info['PoissonRatio_yx'] = nu
  defect.info['OrigHeight']      = orig_height
  
  strain_atoms = fix_edges_defect(defect)
  xpos = defect.positions[:,0]
  ypos = defect.positions[:,1]
  seed = 0.0
  tip  = 0.0
  strain = G_to_strain(initial_G, E, nu, orig_height)
  print 'Applied strain: ', strain*100., '%'
  defect.positions[:,1] += thin_strip_displacement_y(xpos, ypos, strain, seed, tip)
  defect.set_cutoff(3.0)
  defect.calc_connect()
  
  xlen = defect.positions[:,0].max() - defect.positions[:,0].min()
  ylen = defect.positions[:,2].max() - defect.positions[:,2].min()
  print xlen*ylen, orig_height
  A = xlen*ylen
  print 'Applied stress: GPa', E*((strain)/orig_height)/units.GPa

  defect_object             = {}
  defect_object['x']        = x
  defect_object['y']        = y
  defect_object['z']        = z
  defect_object['E']        = E
  defect_object['nu']       = nu
  defect_object['strain_y'] = strain
  defect_object['stress_y'] = E*(strain/orig_height)/units.GPa
  defect_object['T']        = sim_T/units.kB

  defect_json(**defect_object)

  if run_dyn:
    MaxwellBoltzmannDistribution(defect, 2.0*sim_T)
    dynamics = VelocityVerlet(defect, timestep)
    dynamics.attach(pass_print_context(defect, dynamics))
    trajectory = AtomsWriter('{0}_traj.xyz'.format(input_file))
    dynamics.attach(trajectory, print_interval, defect)
    print 'Running Crack Simulation'
    dynamics.run(nsteps)
    print 'Crack Simulation Finished'

  if calc_nye:
    ats = AtomsReader('{0}_traj.xyz'.format(input_file))
    hyb = np.zeros(ats[0].n)
    print len(hyb)
    print ats[0].params
#calc_nye tensor for each configuration in the trajectory:
    print len(ats)
    traj_burg = AtomsWriter('{0}_burg.xyz'.format(input_file))
    for i, at in enumerate(ats[::5]):
      print i
      pp_nye_tensor(at)
      traj_burg.write(at)
