"""
Module containing DistributedForceMixingPotential class.
Calculations are distributed to a number of compute nodes, controlled
by an AtomsServer created by the DistributedForceMixingPotential.
"""
import os
import re
import threading

import numpy as np

from quippy.system import system_timer
from quippy.potential import Potential, ForceMixingPotential
from quippy.clusters import (HYBRID_ACTIVE_MARK, HYBRID_NO_MARK, HYBRID_BUFFER_MARK,
                             construct_hysteretic_region,
                             create_hybrid_weights,
                             create_cluster_simple)
from quippy.farray import farray, fzeros, frange, FortranArray, fenumerate
from quippy.ringstat import distance_map
import quippy.util
#from atomsserver import AtomsServer, AtomsRequestHandler
from matscipy.socketcalc import AtomsServerAsync, AtomsRequestHandler

class DistributedForceMixingPotential(ForceMixingPotential):
    """
    Subclass of ForceMixingPotential which uses separate
    clients to compute MM and QM forces, and communicates with
    them via sockets using an AtomsServer.
    """
    def __init__(self, mm_clients, qm_clients, ip, port=0, rundir=None, bulk_scale=None, mpi_obj=None,
                 callback=None, calculator=None, cutoff_skin=1.0, atoms=None,
                 qm_list=None, fpointer=None, finalise=True, error=None, cluster_args=None, test_mode=False,
                 save_clusters=False, force_restart=False, **kwargs):

        def callback_factory(force_prop_name):
            def callback(at):
                at.add_property('force', getattr(at, force_prop_name),
                                overwrite=True)
            return callback

        if isinstance(mm_clients, Potential):
            self.mm_local = True
            self.pot1 = mm_clients
            self.mm_clients = []
        else:
            self.mm_local = False
            self.mm_clients = mm_clients

        if isinstance(qm_clients, Potential):
            self.qm_pot = qm_clients
            self.qm_clients = []
        else:
            self.qm_clients = qm_clients
        
        clients = self.mm_clients + self.qm_clients
        print 'Clients', clients
        if len(clients) > 0:
            self.server = AtomsServerAsync((ip, port), AtomsRequestHandler, clients, bgq=True)
            self.server_thread = threading.Thread(target=self.server.serve_forever)
            self.server_thread.daemon = True
            self.server_thread.start()

        if rundir is None:
            rundir = os.getcwd()
        self.rundir = rundir

        self.label = 1
        self.cluster_args = {}
        if cluster_args is not None:
            self.cluster_args = cluster_args

        if not self.mm_local:
            self.pot1 = Potential('CallbackPot', callback=callback_factory('mm_force'))
        self.pot2 = Potential('CallbackPot', callback=callback_factory('qm_force'))

        self.test_mode = test_mode
        self.save_clusters = save_clusters
        self.force_restart = force_restart

        ForceMixingPotential.__init__(self, self.pot1, self.pot2, bulk_scale=bulk_scale,
                                      mpi_obj=mpi_obj, cutoff_skin=cutoff_skin,
                                      atoms=atoms, qm_list=qm_list, 
                                      fpointer=fpointer, finalise=finalise,
                                      error=error, **kwargs)

    def calc(self, at, energy=None, force=None, virial=None,
             local_energy=None, local_virial=None,
             args_str=None, error=None, **kwargs):

        clusters = []
        orig_label = self.label
        if not self.mm_local:
            # always submit the MM calc
            if self.test_mode:
                clusters.append(at)
            else:
                self.server.put(at, 0, self.label, force_restart=self.force_restart)
            self.label += 1
        
        do_qm = not self.get('method').startswith('lotf') or self.get('lotf_do_qm')
        
        if do_qm:
            #print 'REGIONS', [ k for k in at.properties.keys() if re.match('hybrid_[0-9]+', k) ]
            n_region = len([ k for k in at.properties.keys() if re.match('hybrid_[0-9]+', k) ])

            if self.get('calc_weights'):
                system_timer('create_hybrid_weights')
                # overall hybrid property is union of all the hybrids
                if not hasattr(at, 'hybrid'):
                    at.add_property('hybrid', 0)
                at.hybrid[:] = 0
                for i in frange(n_region):
                    hybrid = getattr(at, 'hybrid_%d' % i)
                    at.hybrid[hybrid == HYBRID_ACTIVE_MARK] = HYBRID_ACTIVE_MARK
                if not hasattr(at, 'hybrid_mark'):
                    at.add_property('hybrid_mark', at.hybrid)
                at.hybrid_mark[:] = 0
                at.hybrid_mark = at.hybrid
                create_hybrid_weights_args = self.cluster_args.copy()
                create_hybrid_weights_args['buffer_hops'] = 0
                create_hybrid_weights_args['transition_hops'] = 0 # ensure a fast exit
                # overall hybrid -> hybrid_mark, weight_region1
                create_hybrid_weights(at, args_str=quippy.util.args_str(create_hybrid_weights_args))
                system_timer('create_hybrid_weights')
            
            # make clusters and submit to QM clients
            system_timer('make_clusters')
            for i in frange(n_region):
                hybrid_name = 'hybrid_%d' % i
                hybrid_mark_name = 'hybrid_mark_%d' % i
                if self.get('calc_weights'):
                    hybrid = getattr(at, hybrid_name)
                    if not hasattr(at, hybrid_mark_name):
                        at.add_property(hybrid_mark_name, HYBRID_NO_MARK)
                    hybrid_mark = getattr(at, hybrid_mark_name)

                    # set marks to allow previously active atoms to become buffer atoms
                    # create_hybrid_weights will then set the buffer marks
                    hybrid_mark[hybrid_mark == HYBRID_ACTIVE_MARK] = HYBRID_BUFFER_MARK
                    hybrid_mark[hybrid == HYBRID_ACTIVE_MARK] = HYBRID_ACTIVE_MARK

                    print ('region %d, sum(hybrid) %d, sum(hybrid_mark) %d' % 
                            (i, sum(hybrid), sum(hybrid_mark)))

                    create_hybrid_weights_args = self.cluster_args.copy()
                    create_hybrid_weights_args['run_suffix'] = '_%d' % i
                    create_hybrid_weights_args_str = quippy.util.args_str(create_hybrid_weights_args)
                    print 'calling create_hybrid_weights with args_str %s' % create_hybrid_weights_args_str
                    create_hybrid_weights(at, args_str=create_hybrid_weights_args_str)

                cluster_args_str = quippy.util.args_str(self.cluster_args)
                print 'calling create_cluster_simple with args_str %s' % cluster_args_str
                c = create_cluster_simple(at, mark_name=hybrid_mark_name, args_str=cluster_args_str)
                
                client_id = i
                if self.mm_local:
                    client_id -= 1
                if self.save_clusters:
                    c.write(os.path.join(self.rundir,
                                         'cluster.client-%03d.label-%04d.xyz' %
                                         (client_id, self.label)))
                if self.test_mode:
                    clusters.append(c)
                else:
                    self.server.put(c, client_id, self.label, force_restart=self.force_restart)
                self.label += 1
            system_timer('make_clusters')

        # wait for results to be ready
        system_timer('get_results')
        if self.test_mode:
            results = []
            for i, cluster in enumerate(clusters):
                result = cluster.copy()
                result.set_cutoff(self.qm_pot.cutoff())
                result.calc_connect()
                self.qm_pot.calc(result, force=True)
                result.params['label'] = orig_label + i
                results.append(result)
        else:
            results = self.server.get_results()
        system_timer('get_results')

        system_timer('process_results')
        if self.mm_local:
            qm_results = results
        else:
            # process MM results
            mm_results, qm_results = results[:1], results[1:]
            mm_result = mm_results[0]
            at.add_property('mm_force', mm_result.force, overwrite=True)

        if do_qm:
            # process QM results
            at.add_property('qm_force', 0., n_cols=3, overwrite=True)

            # extract forces from each cluster
            for i, r in fenumerate(qm_results):
                #print 'qm forces (orig order?):'
                #print '\n'.join(['%d: pos=%s f=%s' % (j, p, f) for j, p, f in zip(r.index, r.pos, r.force)])

                client_id = i
                if self.mm_local:
                    client_id -= 1
                if self.save_clusters:
                    r.write(os.path.join(self.rundir,
                                         'results.client-%03d.label-%04d.xyz' %
                                         (client_id, r.label)))

                mark_name = 'hybrid_mark_%d' % i
                mask = getattr(r, mark_name) == HYBRID_ACTIVE_MARK
                #print 'Cluster %d: QM forces on atoms %s' % (i, r.index[mask])
                #print r.force[:, mask].T
        # HL if set_fortran is false we need to reduce the index here because
        # the atoms object is expecting python indexing.
                at.qm_force[:, [ind-1 for ind in list(r.index[mask])]] = r.force[:, mask]

        system_timer('process_results')

        # now call the parent calc() to do the force mixing etc.
        force_mixing_args = kwargs.copy()
        force_mixing_args.update(self.cluster_args)
        force_mixing_args['calc_weights'] = False # already done this above
        #print 'calling ForceMixingPotential.calc() with args %r' % force_mixing_args
        ForceMixingPotential.calc(self, at, energy, force, virial, local_energy, 
                                  local_virial, args_str, error, **force_mixing_args)

