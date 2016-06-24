import ase.units
from   quippy.io import AtomsWriter



# For isolated LOTF Defect calculation we will typically want to 
# strain the unit cell and watch how the defect migrates.
# Initial Parameters:
sim_T      = 300.0*units.kB
input_file = 'si111111.xyz'
traj_file  = 's111111_traj.xyz'
timestep   =  1.0*units.fs,   # Timestep (NB: time base units are not fs!)
initial_G  = 4.5*(units.J/units.m**2), # Initial energy flow to crack tip

#Create unit cell with the orientation:
x = [1,1,-2]
y = [-1,1,0]
z = [1,1,1]

screw_slab_unit = BodyCenteredCubic(directions = [x, y, z], size=(1,1,1),
                                    symbol='Fe', pbc=(1,1,1), latticeconstant=2.83)

POT_DIR  = '/Users/lambert/pymodules/imeall/imeall/potentials'
pot_file = os.path.join(POT_DIR, 'Fe_Mendelev.xml')
pot      = Potential('IP EAM_ErcolAd ', param_filename=pot_file)


cij = pot.get_elastic_constants(screw_slab_unit)
print ((cij / units.GPa).round(2))
E   = youngs_modulus(cij, y)
nu  = poisson_ratio(cij, y, x)
print E, nu


#Load atoms and set potential
defect = Atoms(input_file)
defect.set_calculator(pot)

top    = defects.positions[:,1].max()
bottom = defects.positions[:,1].min()
left   = defects.positions[:,0].min()
right  = defects.positions[:,0].max()

strain_atoms = fix_edges(defect)

xpos = crack_slab.positions[:,0]
ypos = crack_slab.positions[:,1]
crack_slab.positions[:,1] += strain

seed = 0.0
tip  = L_x -70.

MaxwellBoltzmannDistribution(defect, 2.0*sim_T)
dynamics = VelocityVerlet(defect, timestep)

dynamics.attach(pass_print_context(defect, dynamics, mm_pot))
dynamics.attach(check_if_cracked_context(strain_atoms,mm_pot), 1, defect)

trajectory = AtomsWriter(traj_file)
dynamics.attach(trajectory, print_interval, atoms)
print 'Running Crack Simulation'
dynamics.run(nsteps)
print 'Crack Simulation Finished'

