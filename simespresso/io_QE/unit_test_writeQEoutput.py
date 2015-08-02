__author__ = 'jeremy'

__author__ = 'jeremy'

import unittest
import espresso_data_file_read
import subprocess
import logging
import sys
from simphony.core.data_container import DataContainer
from simphony.core.cuba import CUBA
from simphony.core.keywords import KEYWORDS
from simphony.cuds.particles import Particle, Particles
import espresso_data_file_write

class OutcomesTest(unittest.TestCase):


    def test_espresso_data_file_write(self):
        print('starting test of data file handler')
        espresso_input_filename = 'test_pw.in'
        SP = DataContainer()
        pc = Particles('quantum_espresso_particles')
        #write parameters for a particular working input file
        SP[CUBA.TORQUE] = 'scf'  #calculation
        SP[CUBA.ZETA_POTENTIAL] = 'from_scratch' #restart_mode
        SP[CUBA.YOUNG_MODULUS] = './'  #pseudo dir
        SP[CUBA.VOLUME_FRACTION] = 'NANO'  #prefix
        SP[CUBA.AMPHIPHILICITY] = '.true.' #tprnfor
        SP[CUBA.NUMBER_OF_TIME_STEPS] = 82000 #max_seconds
        SP[CUBA.DIRECTION] = './' #outdir
###        SP[CUBA.ROLLING_FRICTION] = '.true.' #tprnfor
        SP[CUBA.ROLLING_FRICTION] = 8 #ibrav
        SP[CUBA.ORIGINAL_POSITION] = [0,1,2] #write to celldm(1),(2),(3)
        #write n_atoms  nat
        #write n_atom_types  ntyp
        SP[CUBA.LN_OF_RESTITUTION_COEFFICIENT] = 60.0 #ecutwfc
        SP[CUBA.POISSON_RATIO] = 240 #ecutrho
        SP[CUBA.LATTICE_SPACING] = 'vdw-df-c09' #input_dft
        SP[CUBA.SMOOTHING_LENGTH] = 'local-TF' #mixing_mode
        SP[CUBA.PHASE_INTERACTION_STRENGTH] = 0.8 #mixing_beta
        SP[CUBA.DEBYE_LENGTH] = '1.0d-7' #conv_thr
        SP[CUBA.CHEMICAL_SPECIE] = ['C']
        SP[CUBA.MASS] = [12.0107]
        SP[CUBA.FRICTION_COEFFICIENT] = ['06-C.GGA.fhi.UPF']
        SP[CUBA.PROBABILITY_COEFFICIENT] = 'automatic' #mode
        SP[CUBA.EQUATION_OF_STATE_COEFFICIENT] = [1,4,1,0,0,0] #K_POINTS
        SP[CUBA.KINEMATIC_VISCOSITY] = '(angstrom)'#ATOMIC_POSITIONS

#this actually should go below except there is a'pacticles has no attribute add_particle' problem
        espresso_data_file_write.WriteEspressoInputFile(espresso_input_filename,SP,pc)
        if(0):
            p1 = Particle([1.0,2.0,3.0])
            p1.data[CUBA.CHEMICAL_SPECIE] = 'C'
            pc.add_particle(p1)
            p2 = Particle([2.0,3.0,4.0])
            p2.data[CUBA.CHEMICAL_SPECIE] = 'C'
            pc.add_particle(p2)
            p3 = Particle([3.0,4.0,5.0])
            p3.data[CUBA.CHEMICAL_SPECIE] = 'C'
            pc.add_particle(p3)
            p4 = Particle([4.0,5.0,6.0])
            p4.data[CUBA.CHEMICAL_SPECIE] = 'C'
            pc.add_particle(p4)

            espresso_data_file_write.WriteEspressoInputFile(file_name,SP,pc)

if __name__ == '__main__':
    unittest.main()


