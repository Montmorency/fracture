#fracture

`fracture` is a collection of automated workflow 
scripts for fracture simulations using the
[QUIP](https://libatoms.github.io/QUIP/quippy.html)
package.

The central object is the CrackCell() object defined in crack.py. This
object defines the crack cell geometry, the initial energy flow to the crack
tip, and the type of potential to be used.

Each job folder needs to contain a pickled dictionary object
with the CrackCell parameters. This dictionary can be generated using 
make_dict.py and can be modified by hand.  (TODO: Switch this to json!)

This initializes a crack simulation cell, and outputs
the initial strain conditions, and the surface energy of
the intended cleavage plane, it also generates a file 'crack.xyz'
with the relaxed initial crack_slab.

simulate_crack.py contains the helper functions that are attached to
an LOTFDynamics object i.e. the traj_writer, the qm/mm potentials etc.
The default start file name is 'crack.xyz' default output is 'traj_lotf_2.xyz'.

## Typical Usage
Generate a CrackCell dictionary:
  python make_dict.py

For intergranular fracture run:
  python gb_gen.py

To initiate a seed from the crack_dictionary:
  python crack.py

To begin the dynamical simulation:
  python run_crack.py


