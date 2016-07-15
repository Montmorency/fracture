import os
import sys
import ase.units as units
try:
  from atomsserver import QUIPClient, VaspClient
  from bgqtools import (get_bootable_blocks, boot_blocks, block_corner_iter, 
                        get_hostname_ip, get_cobalt_info, set_unbuffered_stdout)
  BGQAVAIL = True
except:
  print 'No Bluegene tools'
  BGQAVAIL = False

test_mode = len(sys.argv[1:]) > 0 and sys.argv[1] == '-t'
# ******* Start of parameters ***********

input_file     = 'e111112.xyz'     # starting configuration
reference_file = 'ref_slab.xyz' # Reference file for Nye tensor
continuation = True             # If true, restart form last frame of most recent *.traj.xyz file
classical    = False             # If true, do classical MD instead of QM/MM
sim_T        = 300.0*units.kB    # Simulation temperature
rescale_velo = False             # Rescale velocities to 2*sim_T  
timestep     = 2.0*units.fs      # Timestep (NB: time base units are not fs!)
cutoff_skin  = 2.0*units.Ang     # Amount by which potential cutoff is increased
                                 # for neighbour calculations
traj_file = '%s.traj.xyz'        # Trajectory output file
#traj_interval = 10              # Number of time steps between
                                 # writing output frames
param_file   = 'Fe_Mendelev.xml'   # Filename of XML file containing
                                 # potential parameters
mm_init_args = 'IP EAM_ErcolAd'  # Initialisation arguments for
                                 # classical potential

# additional parameters for the QM/MM simulation:
qm_inner_radius   = 3.0*units.Ang # Inner hysteretic radius for QM region
qm_outer_radius   = 5.0*units.Ang # Outer hysteretic radius for QM region
hyst_buffer_inner = 4.0 # Inner hysteretic radius for QM region
hyst_buffer_outer = 6.0 # Outer hysteretic radius for QM region
extrapolate_steps = 5   # Number of steps for predictor-corrector
                        # interpolation and extrapolation
traj_interval = extrapolate_steps
check_force_error = False         # Do QM calc at each step to check pred/corr. err
nsteps = 1000

save_clusters = False  # if True, write each QM cluster to .xyz files
force_restart = True  # if True, force VASP to restart for each QM cluster

traj_index = 1
traj_file = '%d.traj.xyz' % traj_index
while os.path.exists(traj_file):
    traj_index += 1
    traj_file = '%d.traj.xyz' % traj_index

cluster_args = dict(single_cluster=False,
                    cluster_calc_connect=False,
                    cluster_hopping=False,
                    cluster_hopping_nneighb_only=True,
                    cluster_periodic_z = True, # match Gamma vs. kpts
                    cluster_vacuum     = 5.0,
                    hysteretic_buffer=True,
                    hysteretic_buffer_inner_radius=hyst_buffer_inner,
                    hysteretic_buffer_outer_radius=hyst_buffer_outer,
                    min_images_only=True,
                    terminate=False,
                    force_no_fix_termination_clash=True,
                    randomise_buffer=False)
# VASP arguments
#vasp_args=dict(xc='PBE', amix=0.1, amin=0.1, bmix=0.001, amix_mag=0.1, bmix_mag=0.0001, 
#               kpts=[1, 1, 7], kpar=4, lreal='auto', ibrion=13, nsw=1000000, nelmdl=-15, ispin=2,
#               nelm=100, algo='VeryFast', npar=32, lplane=False, lwave=False, lcharg=False, istart=0,
#               voskown=1, ismear=1, sigma=0.1, isym=0) # possibly try iwavpr=12, should be faster if it works

#HL check if magnetic moments work equally well.
vasp_args=dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001, 
               kpts=[1, 1, 7], kpar=4, lreal='auto', ibrion=13, nsw=1000000, nelmdl=-15, ispin=2,
               nelm=100, algo='VeryFast', npar=32, lplane=False, lwave=False, lcharg=False, istart=0,
               voskown=1, ismear=1, sigma=0.1, isym=0) # possibly try iwavpr=12, should be faster if it works

n_core = 1 # number of quantum regions

# Blue Gene specific parameters
acct = 'SiO2_Fracture'
runtime = 360
queue = 'default'
n_qm_jobs = n_core
njobs = n_qm_jobs

#qm_exe = '/home/fbianchi/vasp5/vasp.5.3.new/vasp.bgq'
qm_exe = '/projects/SiO2_Fracture/iron/vasp.bgq'
#qm_exe = '/home/fbianchi/project/exe/vasp5.O3.cplx.sock'
qm_npj = 512
qm_ppn = 1

nodes = qm_npj*njobs
rundir = os.getcwd()
qm_env = {'OMP_NUM_THREAD': 1}

# ***** Setup clients and server *******
if BGQAVAIL and not test_mode:
    hostname, ip =  get_hostname_ip()
    try:
        partsize, partition, job_name = get_cobalt_info()
        jobname = '%s.%s' % (job_name, os.path.splitext(os.path.basename(sys.argv[0]))[0])
    except KeyError:
        # Not running under cobalt, so let's qsub ourselves
        qsub_args = 'qsub -A %s -n %d -t %d -q %s --mode script --disable_preboot %s' % \
            (acct, nodes, runtime, queue, ' '.join(sys.argv))
        print qsub_args
        os.system(qsub_args)
        sys.exit(1)

    blocks = get_bootable_blocks(partition, nodes)
    print('Available blocks: %s' % blocks)
    boot_blocks(blocks)

    qm_subblocks = [(i, bcs) for (i, bcs) in enumerate(block_corner_iter(blocks, qm_npj)) ]
    print 'qm_subblocks', qm_subblocks

    qm_clients = [VaspClient(client_id, qm_exe, qm_env, qm_npj, qm_ppn, block, corner, shape, jobname, 
                             **vasp_args) for client_id, (block, corner, shape) in qm_subblocks ]

else:
    qm_clients = []
    hostname, ip = '<dummy>', 0

print 'FEN hostname: %s' % hostname
print 'FEN server IP: %s' % ip
print 'QM nodes per job: %r' % qm_npj
print 'QM MPI tasks per node: %r' % qm_ppn
print 'Number of QM jobs: %d' % n_qm_jobs
print 'Total number of sub-block jobs: %d' % njobs
print 'Total number of nodes: %d' % nodes

# ******* End of parameters *************
