import os
import json
import argparse

import ase.units as units
from   ase.optimize import FIRE
from   ase.md.verlet import VelocityVerlet
from   ase.constraints import FixAtoms
from   ase.lattice.cubic import Diamond, BodyCenteredCubic
from   ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from   quippy            import Atoms, Potential, calc_nye_tensor, set_fortran_indexing
from   quippy.io         import AtomsReader, AtomsWriter
from   quippy.lotf       import LOTFDynamics, update_hysteretic_qm_region
from   quippy.crack      import thin_strip_displacement_y, G_to_strain, get_strain, strain_to_G
from   quippy.system     import verbosity_push, PRINT_VERBOSE 
from   quippy.potential  import ForceMixingPotential, Potential
from   quippy.elasticity import youngs_modulus, poisson_ratio, rayleigh_wave_speed, AtomResolvedStressField

from simulate_crack import pass_print_context
from tb_pot         import TightBindingPot

import numpy as np

import params
from run_qmmm_111screw import set_quantum, update_qm_region

set_fortran_indexing=False

def traj_writer(dynamics):
  if params.extrapolate_steps == 1 or dynamics.state == LOTFDynamics.Interpolation:
    trajectory.write(dynamics.atoms)

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

def pp_nye_tensor(at, rr=10.0, dis_type='edge'):
# Post Processing routine to append nye tensor information
# to atoms object.
# Load reference slab and calculate connectivity
    ref_slab   = Atoms('./ref_slab.xyz')
    at.set_cutoff(3.0)
    at.calc_connect()
    ref_slab.set_cutoff(3.0)
    ref_slab.calc_connect()
    core       = np.zeros(3)
# To save time it is possible to only
# calculate the nye tensor for a 'core' region
# determined by cutting out a cluster:
    try:
      core   = at.info['core']
    except:
      print 'No Core info'
      core   = np.array([98.0, 98.0, 1.49]) 
    print '\t Core: ', core
    fixed_mask = (np.sqrt((at.positions[:,0]-core[0])**2 + (at.positions[:,1]-core[1])**2) < rr)
    at.add_property('edgex', 0.0)
    at.add_property('edgey', 0.0)
    at.add_property('screw', 0.0)
    cl = at.select(mask=fixed_mask, orig_index=True) 
    print '\t Size of cluster: ', len(cl)
# Append screw and edge information
# Probably should just pass an array at this stage?
    alpha     = calc_nye_tensor(cl, ref_slab, 3, 3, cl.n)
    cl.edgex  = alpha[2,0,:]
    cl.edgey  = alpha[2,1,:]
    cl.screw  = alpha[2,2,:]
# Update quantum region according to the
# position of the dislocation core:
    if dis_type  == 'screw':
      defect_pos = cl.screw
    elif dis_type == 'edge':
      defect_pos = cl.edgex
    total_def    = 0.
    c            = np.array([0.,0.,0.])
    for i in range(len(cl)):
        total_def = total_def + defect_pos[i]
        c[0] = c[0] + cl.positions[i,0]*defect_pos[i]
        c[1] = c[1] + cl.positions[i,1]*defect_pos[i]
    c[0] = c[0]/total_def
    c[1] = c[1]/total_def
    c[2] = cl.lattice[2,2]/2.
    core[:] = c.copy()
# Now thread atoms back in looks like select will have indices in the fortran convention.
    for index, screw, edgex, edgey in zip(cl.orig_index, cl.screw, cl.edgex, cl.edgey):
      at.screw[index-1] = screw
      at.edgex[index-1] = edgex
      at.edgey[index-1] = edgey
    return core
# Initial Parameters:
input_file = 's111111.xyz'
traj_file  = 's111111_traj.xyz'
timestep   =  2.0*units.fs #Timestep (NB: time base units are not fs!)
print_interval = 100
initial_G      = 0.0*(units.J/units.m**2) #Initial energy flow to crack tip or strain energy of slab
nsteps         = 500  # Total number of timesteps to run for

#Create unit cell with the orientation:

x = [1,1,-2]
y = [-1,1,0]
z = [1,1,1]
screw_slab_unit = BodyCenteredCubic(directions = [x, y, z], size=(1,1,1),
                                    symbol='Fe', pbc=(1,1,1), latticeconstant=2.83)
if __name__=='__main__':
#Parse Arguments:
  parser = argparse.ArgumentParser()
  parser.add_argument("-rd",  "--run_dyn", action='store_true')
  parser.add_argument("-cn",  "--calc_nye", action='store_true')
  parser.add_argument("-inp", "--input_file", required=True)
  parser.add_argument("-dt",   "--dis_type", required=True)
  parser.add_argument("-pt",  "--pot_type", help='Potential associated with dynamics: TB or EAM', default='EAM')
  parser.add_argument("-st",  "--sim_T", help='Simulation Temperature in Kelvin. Default is 300 K.', type=float, default=300.0)

  args        = parser.parse_args()
  run_dyn     = args.run_dyn
  calc_nye    = args.calc_nye
  input_file  = args.input_file
  dis_type    = args.dis_type
  pot_type    = args.pot_type
  sim_T       = args.sim_T
#set temperature
  sim_T       = sim_T*units.kB
  calc_elastic_constants = False
  dyn_type = 'eam'

  print '\tDislocation Type: ',   dis_type, ' dislocation.'
  print '\tCalculate Dynamics: ', run_dyn
  print '\tCalculate Nye Tensor: ', calc_nye
  print '\tInput File: ', input_file
  print '\tPotential: ',  pot_type

  if pot_type == 'EAM':
    POT_DIR  = '/Users/lambert/pymodules/imeall/imeall/potentials'
    eam_pot  = os.path.join(POT_DIR, 'Fe_Mendelev.xml')
    r_scale  = 1.00894848312
    pot      = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
  elif pot_type == 'TB':
    tb = TightBindingPot(alat= 1.00, nk=12)
    screw_slab_unit = Atoms(screw_slab_unit)
    center_cell = np.diag(screw_slab_unit.get_cell())/2.
    screw_slab_unit.add_atoms(center_cell,1)
    #screw_slab_unit.add_atoms(center_cell*0.76, 6)

    screw_slab_unit.set_atoms(screw_slab_unit.Z)
    print len(screw_slab_unit)
    tb.write_control_file('ctrl.fe', screw_slab_unit)
    pot = Potential('IP LMTO_TBE', param_str="""
    <params>
      <LMTO_TBE_params n_types="2" control_file="ctrl.fe">
        <per_type_data type="1" atomic_num="26"/>
        <per_type_data type="2" atomic_num="1"/>
      </LMTO_TBE_params>
    </params>""")

#If there is an initial strain energy increment that here:
  pot.print_()
  screw_slab_unit.write('screw_unitcell.xyz')
  screw_slab_unit.set_calculator(pot)


#Elastic constants are disabled until we have tight binding operational.
  if calc_elastic_constants and dyn_type != 'LOTF':
    cij = pot.get_elastic_constants(screw_slab_unit)
    print 'ELASTIC CONSTANTS'
    print ((cij / units.GPa).round(2))
    E   = youngs_modulus(cij, y)
    nu  = poisson_ratio(cij, y, x)
    print 'Youngs Modulus: ', E/units.GPa,'Poisson Ratio: ', nu
    print 'Effective elastic modulus E: ', E/(1.-nu**2)

  print 'Loading Structure File:', '{0}.xyz'.format(input_file)
  defect = Atoms('{0}.xyz'.format(input_file))

  top    = defect.positions[:,1].max()
  bottom = defect.positions[:,1].min()
  left   = defect.positions[:,0].min()
  right  = defect.positions[:,0].max()
  orig_height = (defect.positions[:,1].max()-defect.positions[:,1].min())

#Attaching Properties to atoms object
  if calc_elastic_constants:
    defect.info['YoungsModulus']   = E
    defect.info['PoissonRatio_yx'] = nu

  defect.info['OrigHeight']        = orig_height
  strain_atoms = fix_edges_defect(defect)
  xpos = defect.positions[:,0]
  ypos = defect.positions[:,1]
  seed = 0.0
  tip  = 0.0

  if False:
    strain = G_to_strain(initial_G, E, nu, orig_height)
    print 'Applied strain: ', strain*100., '%'
    defect.info['strain'] = strain
    defect.info['G']      = initial_G
    defect.positions[:,1] += thin_strip_displacement_y(xpos, ypos, strain, seed, tip)

  defect.set_cutoff(3.0)
  defect.calc_connect()
  
  xlen = defect.positions[:,0].max() - defect.positions[:,0].min()
  ylen = defect.positions[:,2].max() - defect.positions[:,2].min()

  defect_object             = {}
  defect_object['x']        = x
  defect_object['y']        = y
  defect_object['z']        = z
  defect_object['T']        = sim_T/units.kB
  defect_json(**defect_object)

  if run_dyn:
  # Load atoms and set potential
  # If tight-binding need to re-initialize potential
  # with new defect atoms object.
    if pot_type=='TB':
      tb.write_control_file('ctrl.fe', defect)
      pot = Potential('IP LMTO_TBE', param_str="""
      <params>
        <LMTO_TBE_params n_types="3" control_file="ctrl.fe">
          <per_type_data type="1" atomic_num="26"/>
          <per_type_data type="2" atomic_num="6"/>
          <per_type_data type="3" atomic_num="1"/>
        </LMTO_TBE_params>
      </params>""")

      screw_slab_unit.rattle(0.01)
      opt = FIRE(screw_slab_unit)
      opt.run(fmax=0.01)
      screw_slab_unit.write('screw_unitcell_H.xyz')

      print 'Initializing EAM Potential'
      POT_DIR = '/Users/lambert/pymodules/imeall/imeall/potentials'
      eam_pot = os.path.join(POT_DIR, 'Fe_Mendelev.xml')
      r_scale = 1.00894848312
      mm_pot  = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
  
      print 'Initializing TightBinding Potential'
      tb              = TightBindingPot(alat= 2.87, nk=1)
    #screw_slab_unit = Atoms('clusters.xyz')
      tb.write_control_file('ctrl.fe', screw_slab_unit)
      qm_pot = Potential('IP LMTO_TBE', param_str="""
      <params>
        <LMTO_TBE_params n_types="3" control_file="ctrl.fe">
          <per_type_data type="1" atomic_num="26"/>
          <per_type_data type="2" atomic_num="6"/>
          <per_type_data type="3" atomic_num="1"/>
        </LMTO_TBE_params>
      </params>""")
      print 'Initializing LOTF Potential, qm_radii:', params.qm_inner_radius, params.qm_outer_radius
      qmmm_pot = ForceMixingPotential(pot1=mm_pot, pot2=qm_pot, atoms=defect,
                                 qm_args_str='single_cluster cluster_periodic_z carve_cluster '+
                                'terminate=F cluster_hopping=F randomise_buffer=F',
                                 fit_hops=2,
                                 lotf_spring_hops=2,
                                 hysteretic_buffer=True,
                                 hysteretic_buffer_inner_radius=params.hyst_buffer_inner,
                                 hysteretic_buffer_outer_radius=params.hyst_buffer_outer,
                                 cluster_hopping_nneighb_only=False,
                                 min_images_only=True)
      defect.set_calculator(qmmm_pot)
    elif pot_type == 'EAM':
      print 'Initializing EAM Potential'
      POT_DIR = '/Users/lambert/pymodules/imeall/imeall/potentials'
      eam_pot = os.path.join(POT_DIR, 'Fe_Mendelev.xml')
      r_scale = 1.00894848312
      pot  = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
      defect.set_calculator(pot)
    else:
      print 'No potential chosen', 1/0

    print 'Finding initial dislocation core positions...'
    try:
      defect.params['core']
    except KeyError:
      defect.params['core'] = np.array([98.0, 98.0, 1.49])

    defect  = set_quantum(defect, params.n_core)
    MaxwellBoltzmannDistribution(defect, 2.0*sim_T)
    if dyn_type =='eam':
       dynamics = VelocityVerlet(defect, timestep)
       dynamics.attach(pass_print_context(defect, dynamics))
    elif dyn_type =='LOTF':
       defect.info['core']= np.array([98.0, 98.0, 1.49])
       print 'Initializing LOTFDynamics'
       verbosity_push(PRINT_VERBOSE)
       dynamics = LOTFDynamics(defect, timestep,
                               params.extrapolate_steps,
                               check_force_error=False)
       dynamics.set_qm_update_func(update_qm_region)
       dynamics.attach(pass_print_context(defect, dynamics))
       dynamics.attach(traj_writer, print_interval, defect)
    else:
      print 'No dyn_type chosen', 1/0
    
    trajectory = AtomsWriter('{0}.traj.xyz'.format(input_file))
    print 'Running Crack Simulation'
    dynamics.run(nsteps)
#Write cooked i.e. thermalized ceel to file.
    defect.set_cutoff(3.0)
    defect.calc_connect()
    new_core = pp_nye_tensor(defect, dis_type=dis_type)
    defect.info['core']= new_core
    defect.write('{0}_therm.xyz'.format(input_file))
    print 'Crack Simulation Finished'

  if calc_nye:
#ats = AtomsReader('{0}_traj.xyz'.format(input_file))
    ats = AtomsReader('{0}.xyz'.format(input_file))
    hyb = np.zeros(ats[0].n)
#calc_nye tensor for each configuration in the trajectory:
    print len(hyb)
    print ats[0].params
    print len(ats)
    traj_burg = AtomsWriter('{0}_burg.xyz'.format(input_file))
    for i, at in enumerate(ats[::1]):
      del at.properties['edgex']
      del at.properties['edgey']
      del at.properties['screw']
      if i > 0:
        at.info['core'] = new_core
      new_core = pp_nye_tensor(at, dis_type=dis_type)
      print i, new_core
      traj_burg.write(at)
