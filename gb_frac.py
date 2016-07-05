import json
import numpy as np

import ase.units as units
from   ase.optimize import FIRE

from quippy import Atoms, Potential
from quippy import set_fortran_indexing
from quippy.io           import AtomsWriter, AtomsReader, write
from quippy import frange, farray, fzeros

from ase.lattice.cubic   import BodyCenteredCubic
from ase.lattice.surface import surface, bcc111,bcc110
from ase.utils.geometry  import get_duplicate_atoms

set_fortran_indexing(False)


def del_atoms(x=None):
  rcut = 2.0
  #x = Atoms('crack.xyz')
  if x == None:
    x = Atoms('1109337334_frac.xyz')
  else:
    pass
  x.set_cutoff(3.0)
  x.calc_connect()
  x.calc_dists()
  rem=[]
  r = farray(0.0)
  u = fzeros(3)
  print len(x)
  for i in frange(x.n):
    for n in frange(x.n_neighbours(i)):
      j = x.neighbour(i, n, distance=3.0, diff=u)
      if x.distance_min_image(i, j) < rcut and j!=i:
        rem.append(sorted([j,i]))
    if i%10000==0: print i
  rem = list(set([a[0] for a in rem]))
  if len(rem) > 0:
    print rem
    x.remove_atoms(rem)
  else:
    print 'No duplicate atoms in list.'
  x.write('crack_nodup.xyz')
  return x

class GBFracture(object):
  def build_surface(self, bp=[-1,1,12], v=[1,1,0]):
    '''
    Build Surface unit cell
    '''
    bpxv = [(bp[1]*v[2]-v[1]*bp[2]), (bp[2]*v[0]-bp[0]*v[2]), (bp[0]*v[1]- v[0]*bp[1])]
    surf_cell = BodyCenteredCubic(directions = [v, bpxv, bp],
                                  size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                                  latticeconstant = 2.85)
    n = 2
    while(surf_cell.get_cell()[2,2]< 20.0 ):
      surf_cell = BodyCenteredCubic(directions = [v, bpxv, bp],
                    size = (1,1,n), symbol='Fe', pbc=(1,1,1),
                    latticeconstant = 2.85)
      n += 1
    surf_cell.center(vacuum=20.0, axis=2)
    return surf_cell

  def build_tilt_sym_frac(self, bp=[-1,1,12], v=[1,1,0], c_space=None):
    bpxv = [(bp[1]*v[2]-v[1]*bp[2]),(bp[2]*v[0]-bp[0]*v[2]),(bp[0]*v[1]- v[0]*bp[1])]
    grain_a = BodyCenteredCubic(directions = [bpxv, bp, v],
                              size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                              latticeconstant = 2.85)
    n_grain_unit = len(grain_a)
    n = 2
# Slightly different from slabmaker since in the fracture
# Simulations we want the grainboundary plane to be normal to y.
    while(grain_a.get_cell()[1,1] < 120.0):
      grain_a = BodyCenteredCubic(directions = [bpxv, bp, v],
                           size = (1,n,2), symbol='Fe', pbc=(1,1,1),
                           latticeconstant = 2.85)
      n += 1
    print '\t {0} repeats in z direction'.format(n)
    grain_b = grain_a.copy()
    grain_c = grain_a.copy()
    print '\t', '{0} {1} {2}'.format(v[0],v[1],v[2])
    print '\t', '{0} {1} {2}'.format(bpxv[0],bpxv[1],bpxv[2])
    print '\t', '{0} {1} {2}'.format(bp[0], bp[1], bp[2])
    if c_space==None:
      s1 = surface('Fe', (map(int, bp)), n)
      c_space = s1.get_cell()[2,2]/float(n) #-s1.positions[:,2].max()
      s2 = surface('Fe', (map(int, v)), 1)
      x_space = s2.get_cell()[0,0] #-s1.positions[:,2].max()
      s3 = surface('Fe', (map(int, bpxv)), 1)
      y_space = s3.get_cell()[1,1] #-s1.positions[:,2].max()
    print '\t Interplanar spacing: ', x_space.round(2), y_space.round(2), c_space.round(2), 'A'
	# Reflect grain b in z-axis (across mirror plane):
    print grain_a.get_cell()[1,1]-grain_a.positions[:,1].max()
    grain_b.positions[:,1]  = -1.0*grain_b.positions[:,1]
    grain_c.extend(grain_b)
    grain_c.set_cell([grain_c.get_cell()[0,0], 2*grain_c.get_cell()[1,1], grain_c.get_cell()[2,2]])
    grain_c.positions[:,1] += abs(grain_c.positions[:,1].min())
    pos = [grain.position for grain in grain_c]
    pos = sorted(pos, key= lambda x: x[2])
    dups = get_duplicate_atoms(grain_c)
 #  now center fracture cell with plenty of vacuum
    grain_c.center(vacuum=10.0,axis=1)
    return grain_c

if __name__=='__main__':
  sym_tilt_110 = [[np.pi*(93.37/180.), np.array([-3., 3., 4.])]]
  gbid = '1109337334'
  or_axis = [1,1,0]
  bp      = [-1,1,12]

##########################################
##First calculate surface energetics:   ##
##########################################
  bulk = BodyCenteredCubic(directions = [[1,0,0],[0,1,0], [0,0,1]],
                           size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                           latticeconstant = 2.85)

  eam_pot = './Fe_Mendelev.xml'
  gb_frac = GBFracture()
  surf_cell = gb_frac.build_surface(bp = sym_tilt_110[0][1])
  pot     = Potential('IP EAM_ErcolAd', param_filename=eam_pot)
  bulk.set_calculator(pot)
  ener_per_atom = bulk.get_potential_energy()/len(bulk)
  surf_cell.set_calculator(pot)
  surf_ener = surf_cell.get_potential_energy()
  cell  = surf_cell.get_cell()
  A     = cell[0][0]*cell[1][1]
  gamma = (surf_ener- len(surf_cell)*ener_per_atom)/A

  print '2*gamma ev/A2', gamma
  print '2*gamma J/m2',  gamma/(units.J/units.m**2)
  j_dict = {'or_axis':or_axis, 'bp':bp, 'gamma':gamma}
  with open('gbfrac.json','w') as f:
    json.dump(j_dict, f)
  out = AtomsWriter('{0}'.format('{0}_surf.xyz'.format(gbid)))
  out.write(Atoms(surf_cell))
  out.close()
  frac_cell = gb_frac.build_tilt_sym_frac()

#Unit cell for grain boundary fracture cell:
  print frac_cell.get_cell().round(2)

  frac_cell = Atoms(frac_cell)
  frac_cell = del_atoms(frac_cell)

#Relax grainboundary crack cell unit cell:
  pot      = Potential('IP EAM_ErcolAd', param_filename='Fe_Mendelev.xml')
  frac_cell.set_calculator(pot)
  slab_opt          = FIRE(frac_cell)
  slab_opt.run(fmax = (0.02*units.eV/units.Ang))

#Print frac_cell to file:
  out = AtomsWriter('{0}'.format('frac_cell.xyz'.format(gbid)))
  out.write(Atoms(frac_cell))
  out.close()
