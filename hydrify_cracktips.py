import os
import glob 
import shutil
from quippy import Atoms

jobs = glob.glob('T*H1')

# Create new directories which mimic the original
# pristine crack cell but then add hydrogen at the
# crack tip.
for job in jobs
  os.mkdir(job+'H1')
  input  = os.path.join(job, 'crack.xyz')
  pickle = os.path.join(job, 'crack_info.pckl')
  shutil.copy(input, job+'H1')
  shutil.copy(pickle, job+'H1')

for job in jobs:
  os.chdir(job)
  print job
  cr = Atoms('crack.xyz')
  cr.add_atoms([109.17,-1.11,4.275],1)
  cr.set_atoms(cr.numbers)
  cr.write('crack.xyz')
  os.chdir('../')

