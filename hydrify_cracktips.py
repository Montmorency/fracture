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
    cl         = ats.select(fixed_mask, orig_index=True)
    self.write_cluster(cl.copy(), name = 'cluster.xyz')
    1/0
#  Then compute voronoi diagram and add hydrogens around there:
    vor        = spatial.Voronoi(cl.positions, furthest_site=False)
    delaunay   = spatial.Delaunay(cl.positions, furthest_site=False)
    f = open('furthest_site.dat', 'w')
    g = open('fe_sites.dat', 'w')
    #for i in range(n_h):
    for vert in vor.vertices:
      ats.add_atoms(vert, 1)
    for pos in vor.vertices:
      print >> f, pos
    fixed_mask = (np.sqrt(map(sum, map(np.square, ats.positions[:,0:3]-ats.params['CrackPos'][0:3]))) <= rr)
    cl         = ats.select(fixed_mask)
    self.write_cluster(cl.copy(), name = 'clusterH.xyz')
    for pos in cl.positions:
      print >> g, pos
    f.close()
    g.close()

  def append_if_thresh(h_pos):
    pared_h = h_pos[0]
    for h in h_pos[1:]:
      if all([np.linalg.norm(h-h_unique) > 1.6 for h_unique in pared_h]):
        pared_h =  np.vstack((pared_h, h))
    return pared_h

  def hydrogenate_gb(self, gb, d=2.0, n_H=8):
    """
    Given a grain boundary, find a bulk plane to dump hydrogen in,
    and a platelet of hydrogen parallel to the grain boundary. Routine
    returns positions of a suitable "plane" +/- 2A.
    """
    z_bulk     = gb.lattice[2,2]/2.0
# Select the bulk plane:
    fixed_mask = (np.sqrt(np.square(gb.positions[:,2]-z_bulk)) <= d)
    cl         = gb.select(fixed_mask, orig_index=True)
    cl.write('plane.xyz')
    delaunay   = spatial.Delaunay(cl.positions, furthest_site=False)
#http://stackoverflow.com/questions/10650645/python-calculate-voronoi-tesselation-from-scipys-delaunay-triangulation-in-3d?rq=1
    p = tri.points[tri.vertices]
    A = p[:,0,:].T
    B = p[:,1,:].T
    C = p[:,2,:].T
    a = A - C
    b = B - C
    def dot2(u, v):
      return u[0]*v[0] + u[1]*v[1]

    def cross2(u, v, w):
      """u x (v x w)"""
      return dot2(u, w)*v - dot2(u, v)*w

    def ncross2(u, v):
      """|| u x v ||^2"""
      return sq2(u)*sq2(v) - dot2(u, v)**2

    def sq2(u):
      return dot2(u, u)
#answer then goes on to grab the Voronoi edges but we just
#need the circumcenters
#We want to append hydrogens at unique sites in the circumcenters.
#This list has the circumcenters for lots of different simplices
    cc = cross2(sq2(a) * b - sq2(b) * a, a, b) / (2*ncross2(a, b)) + C
    plane_cc = cc.T
    h_pos=filter(lambda x: np.sqrt(np.square(x[2]-z_bulk)) <= d, plane_cc)
#so we removed them
    h_pos = append_if_thresh(h_pos)
    for h in h_pos:
      cl.add_atoms(h,1)
    cl.write('hydrogenated_grain')
    return

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

