import os
import glob 
import shutil
import argparse 
import numpy as np
import scipy.spatial as spatial

from   quippy import Atoms, set_fortran_indexing

set_fortran_indexing(False)

class Hydrify(object):
  """"
  class:Hydrify contains methods for adding hydrogens at specific positions in a 
  simulation cell.
  """ 
  def hydrify_multi_dirs(self, jobs, crack_pos):
    """
    Create new directories which mimic the original
    pristine crack cell but then add hydrogen at the
    crack tip. 
    """
    for job in jobs:
      os.mkdir(job+'H1')
      input  = os.path.join(job, 'crack.xyz')
      pickle = os.path.join(job, 'crack_info.pckl')
      shutil.copy(input, job+'H1')
      shutil.copy(pickle, job+'H1')
    for job in jobs:
      os.chdir(job)
      print job
      cr = Atoms('crack.xyz')
      self.hydrify_single_job(cr)
      cr.add_atoms([109.17,-1.11,4.275],1)
      cr.set_atoms(cr.numbers)
      cr.write('crack.xyz')
      os.chdir('../')

  def hydrify_single_job(self, ats, positions=[]):
    """
    Add hydrogens at positions in list.
    """
    for h_pos in positions:
      ats.add_atoms([109.17,-1.11,4.275],1)
      ats.set_atoms(cr.numbers)

  def write_cluster(self, cl,name='cluster.xyz'):
    cl.center(vacuum=2.0)
    cl.write(name)

  def add_h_smart(self, ats, position=np.array([0.0,0.0,0.0]),rr=10.0,n_h=1):
    """
    Given a position try to add the hydrogen in a "nice"
    place. rr is radius of clust
    """
# First take a small selection of atoms around the crack tip
    ats.calc_connect(2.5)
    h_pos = np.array([118.032, 2.55, 2.78])
    fixed_mask = (np.sqrt(map(sum, map(np.square, ats.positions[:,0:3]-h_pos[0:3]))) <= rr)
    #fixed_mask = (np.sqrt(map(sum, map(np.square, ats.positions[:,0:3]-ats.params['CrackPos'][0:3]))) <= rr)
    cl         = ats.select(fixed_mask, orig_index=True)
    self.write_cluster(cl.copy(), name = 'cluster.xyz')
    1/0
#  Then compute voronoi diagram and add hydrogens around there:
    vor        = spatial.Voronoi(cl.positions, furthest_site=False)
    delaunay   = spatial.Delaunay(cl.positions, furthest_site=False)
    f = open('furthest_site.dat', 'w')
    g = open('fe_sites.dat', 'w')
    #for i in range(n_h):
    for i in range(len(vor.vertices)):
      ats.add_atoms(vor.vertices[i], 1)
    for pos in vor.vertices:
      print >> f, pos
    fixed_mask = (np.sqrt(map(sum, map(np.square, ats.positions[:,0:3]-ats.params['CrackPos'][0:3]))) <= rr)
    cl         = ats.select(fixed_mask)
    self.write_cluster(cl.copy(), name = 'clusterH.xyz')
    for pos in cl.positions:
      print >> g, pos
    f.close()
    g.close()

if __name__=='__main__':
  parser  = argparse.ArgumentParser()
  parser.add_argument('-p','--pattern', default="1.traj.xyz")
  args    = parser.parse_args()
  pattern = args.pattern

  jobs = glob.glob(pattern)
  print jobs
  hydrify =  Hydrify()
  for job in jobs:
    ats = Atoms(job)
    print job
    print len(ats)
    print ats.params['CrackPos']
    hydrify.add_h_smart(ats, position=[ats.params['CrackPos']])
    #hydrify.add_h_smart(ats, position=[h_pos])
    ats.write('thermalizedHtest.xyz')

