Automated workflow for fracture simulation. To be  used in conjunction with 
automated job generation and submission in pwsgwpy.

1) crack.py contains CrackCell() object.
	 Each job folder will contain a pickled dictionary object 
	 with the Crack cell parameters. This dictionary can be modified by hand
	 or generated automatically.

	 This initializes a crack simulation cell, and outputs 
	 the initial strain conditions, and the surface energy of
   the intended cleavage plane.

	 Also generates a file 'crack.xyz' with a relaxed initial crack_slab.

2) simulate_crack.py contains the helper functions that are attached to
	 an LOTFDynamics object i.e. traj_writer, set the qmmm potentials etc.

3) In the submission script we execute:
	 		python run_crack.py

	 Options can pass the script the --restart parameter if we want the
	 calculation to restart from a particular trajectory file, specify the
	 name of the output file with:
				python run_crack.py --restart 'file' --output 'outputfile'

	The default start file name is 'crack.xyz' default output is 'traj.xyz'





