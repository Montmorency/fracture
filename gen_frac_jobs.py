import sys
import os
import shutil
import itertools
import pickle
import ase.units as units

from fracture import crack as frac

# Open a new directory for each
# configuration, input crack_dict for that 
# configuration.	
G               = [4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
cleavage_planes = [(1,1,1)]
pbs = open('calc.pbs', 'r').read()
# Crack Geometry 1
crack_direction = (-2, 1, 1)    # Miller index of x-axis
cleavage_plane  = (1, 1, 1)     # Miller index of y-axis
crack_front     = (0, 1, -1)    # Miller index of z-axis
# System 2. (110)[001]
#crack_direction = (1,-1,0)
#cleavage_plane = (1,1,0)
#crack_front = (0,0,1)
# System 3. (110)[1-10]
#crack_direction = (0,0,-1)
#cleavage_plane  = (1,1,0)
#crack_front     = (1,-1,0)

for cleave, g in itertools.product(*[cleavage_planes, G]):
	cleave_string = '{0}{1}{2}'.format(cleave[0], cleave[1], cleave[2])
	g_string      = str(g).replace('.','')
	dir_string    = cleave_string+'G'+g_string
	print dir_string
#Build the directory
	try:
		os.mkdir(dir_string)
	except OSError as exc:
		print exc
		pass
#Dump a pbs with appropriate name into the directory
	target_dir = './'+dir_string
	print >> open(target_dir+'/calc.pbs','w'), pbs.format(name='Si'+dir_string)
#Dump a crack dict with the appropriate parameters into the dictionary
	crack_dict = frac.init_dict(initial_G= g*(units.J/units.m**2), crack_direction=crack_direction,
															crack_front=crack_front, cleavage_plane=cleave, symbol='Si')
	f = open('./'+dir_string+'/crack_info.pckl','w')
	pickle.dump(crack_dict, f)
	f.close()
#Copy the parameters file 
	shutil.copy('params.xml', target_dir)
	print 'qsub', target_dir+'/calc.pbs'

