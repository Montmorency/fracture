import os
import glob 
import shutil
import argparse 
import numpy as np
import scipy.spatial as spatial

from   quippy import Atoms

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

  def add_h_smart(self, ats, position=np.array([0.0,0.0,0.0]),rr=7.0,n_h=1):
    """
    Given a position try to add the hydrogen in a "nice"
    place. rr is radius of clust
    """
# First take a small selection of atoms around the crack tip
    ats.calc_connect(2.5)
    fixed_mask = (np.sqrt(sum(map(np.square,(ats.positions[:,0:3]-ats.params['CrackPos'][:])))) < rr)
    cl         = ats.select(fixed_mask)
#  Then compute voronoi diagram and add hydrogens around there:
    vor        = spatial.Voronoi(cl.positions, furthest_site=True)
    for i in range(n_h):
      ats.add_atoms(vor.vertices[i], 1)


if __name__=='__main__':
  parser  = argparse.ArgumentParser()
  parser.add_argument('-p','--pattern', default="nopattern")
  args    = parser.parse_args()
  pattern = args.pattern

  jobs = glob.glob(pattern+"*xyz")
  hydrify =  Hydrify()

  for job in jobs:
    ats = Atoms(job)
    print job
    print len(ats)
    print ats.params['CrackPos']
    hydrify.add_h_smart(ats, position=[ats.params['CrackPos']])
    ats.write('thermalizedH.xyz')

