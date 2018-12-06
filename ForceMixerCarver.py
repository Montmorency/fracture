import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.utils.geometry import find_mic
from ase.lattice.cubic import FaceCenteredCubic,BodyCenteredCubic

# just allows for pbc / free
class ForceMixingCarvingCalculator:
    def __init__(self, atoms, qm_list, mm_calc, qm_calc, pbc_type=[False,False,False],\
               vacuum=5., alpha=1., beta=1., buffer_width=3.):
        self.name = 'ForceMixingCarvingCalculator'
        self.mm_calc = mm_calc
        self.qm_calc = qm_calc
        self.alpha = alpha
        self.beta = beta
        self.vacuum = vacuum
        self.buffer = buffer_width
        self.qm_region = qm_list

        self.dft_atoms = None
        self.cell = None

        self.qm_positions = atoms.get_positions()[self.qm_region]
        self.qm_shift = np.zeros(3)
        self.system = None

        # Better way of doing this?
        if isinstance(pbc_type, bool):
            self.pbc=[pbc_type,pbc_type,pbc_type]
        elif isinstance(pbc_type,list):
            self.pbc=pbc_type

        mm_pos = atoms.get_positions()
        qm_pos = mm_pos[self.qm_region]
        qm_com = qm_pos.mean(axis=0) # just to avoid pbc errors..
        qm_rad = np.ceil(np.linalg.norm((qm_pos-qm_com)[:,~np.r_[self.pbc]],axis=1).max())


        self.cell = np.identity(3) * (self.vacuum + qm_rad + self.buffer)*2.
        for i in range(3):
            if self.pbc[i]:
                self.cell[:,i] = atoms.get_cell()[:,i]
        # only count distance along non-pbc directions for buffer
        self.dft_atoms = np.where(np.linalg.norm((mm_pos*self.alpha - \
                    qm_com*self.alpha)[:,~np.r_[self.pbc]],axis=1) <= qm_rad + \
                    self.buffer+1.,True,False)

        for ati in np.argwhere(self.dft_atoms).flatten():
            # Minimum distance from cluster <= buffer to be counted
            dist = np.linalg.norm((mm_pos[ati]-qm_pos)[:,~np.r_[self.pbc]],axis=1).min()
            if dist > self.buffer:
                self.dft_atoms[ati] = False

        self.qm_shift = .5*self.cell.diagonal() -\
                        mm_pos[self.dft_atoms].mean(axis=0)*self.alpha


    def get_forces(self,atoms):
        # Scale MD forces to match QM elasticity
        tot_force = self.mm_calc.get_forces(atoms)/self.beta

        # make system object and center atoms in box
        qm_system = atoms.copy()
        del qm_system.constraints
        qm_force = np.zeros(tot_force.shape)
        qm_system = qm_system[self.dft_atoms]

        qm_positions = qm_system.get_positions()
        scaled_qm_positions = qm_positions*self.alpha + self.qm_shift

        qm_system.set_positions(scaled_qm_positions.copy())
        qm_system.set_cell(self.cell.copy())

        dft_atom_count = np.sum(self.dft_atoms)
        tmp_qmf = self.qm_calc.get_forces(qm_system)
        qm_force[self.dft_atoms] = tmp_qmf

        # Reindexing check - will occur when we have multiple species
        if np.linalg.norm(scaled_qm_positions-qm_system.get_positions()) > 1e-9:
            #print "ATOMS REINDEXING IN QM CALC!"
            remap = np.zeros(dft_atom_count).astype(int)
            for i in range(len(remap)):
                remap[i] = find_mic(np.linalg.norm(qm_system.get_positions()-scaled_qm_positions[i],self.cell,pbc=self.pbc))
            qm_force[self.dft_atoms] = tmp_qmf[remap]

        tot_force[self.qm_region] = qm_force[self.qm_region]
        return tot_force

    def get_potential_energy(self, atoms, force_consistent=False):
        return 0.0

    def get_stress(self, atoms):
        raise NotImplementedError
