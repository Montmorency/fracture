import os
import sys
import glob 
import pickle
import subprocess
import numpy as np
from   ase import units
from   quippy import Atoms, set_fortran_indexing
from   quippy import set_fortran_indexing
from   hydrify_cracktips import Hydrify

set_fortran_indexing(False)
scratch = os.getcwd()
hydrify = Hydrify()

ats = Atoms('crack.xyz')
#ats.params['CrackPos'] = np.array([330, 0.0, 0.0])
with open('crack_info.pckl','r') as f:
  crack_dict= pickle.load(f)
print 'G: {}, H_d: {}, sim_T {}'.format(crack_dict['initial_G']*(units.m**2/units.J),
                                        crack_dict['H_d'], crack_dict['sim_T']/units.kB)
h_list  = hydrify.hydrogenate_gb(ats, mode='CrackTip', d_H=crack_dict['H_d'][0], tetrahedral=True, crackpos_fix=ats.params['CrackPos'])
for h in h_list:
  ats.add_atoms(h,1)

#ats.wrap()
ats.write('crackH.xyz')
ats = None
ats = Atoms('crackH.xyz')
ats.set_cutoff(2.4)
ats.calc_connect()
ats.calc_dists()
filter_mask = (ats.get_atomic_numbers()==1)
h_atoms     = ats.select(filter_mask, orig_index=True)
rem=[]
u = np.zeros(3)
for i in h_atoms.orig_index:
  print 'hindex', i
  print 'nneighbs', ats.n_neighbours(i)
  for n in range(ats.n_neighbours(i)):
    j = ats.neighbour(i, n+1, distance=2.4, diff=u)
    print 'neighb index', j
    if ats.distance_min_image(i,j) < 1.1 and j!=i:
      rem.append(i)
rem = list(set(rem))
if len(rem) > 0:
  print 'Removing {} H atoms'.format(len(rem))
  ats.remove_atoms(rem)
else:
  print 'No H atoms closer than threshold.'
#Now a little housekeeping. In the vicinity of a cracktip
#Delaunay can go a little haywire. We remove any H that is far too close to an Fe atom
# |h-fe| < 1.1. and we remove the vacuum Hs
ats.write('crackH.xyz')
h_atoms = sum((ats.get_atomic_numbers() ==1))
zlen    = ats.get_cell()[2,2]
with open('h.txt','w') as f:
  print >>f,  'There are {} H atoms per {} A'.format(h_atoms, zlen)
