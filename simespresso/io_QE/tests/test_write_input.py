__author__ = 'jeremy'

__author__ = 'jeremy'

import unittest
import subprocess
import logging
import sys

from simphony.core.data_container import DataContainer
from simphony.core.cuba import CUBA
from simphony.core.keywords import KEYWORDS
from simphony.cuds.particles import Particle, Particles
from simespresso.io_QE import espresso_class
#import simespresso.io_QE.espresso_class
#import io_QE.espresso_class


class OutcomesTest(unittest.TestCase):


    def test_espresso_data_file_write(self):
        wrp = espresso_class.qe_functions()
        print('starting test of qe input file write')
        espresso_input_filename = 'pw_generated.in'
        wrp.SP = DataContainer()
        wrp.pc = Particles('quantum_espresso_particles')
        #write parameters for a particular working input file
        SP=wrp.SP
        pc=wrp.pc
        SP[CUBA.TORQUE] = 'scf'  #calculation
        SP[CUBA.ZETA_POTENTIAL] = 'from_scratch' #restart_mode
        SP[CUBA.YOUNG_MODULUS] = './'  #pseudo dir
        SP[CUBA.VOLUME_FRACTION] = 'NANO'  #prefix
        SP[CUBA.AMPHIPHILICITY] = '.true.' #tprnfor
        SP[CUBA.NUMBER_OF_TIME_STEPS] = 82000 #max_seconds
        SP[CUBA.DIRECTION] = './' #outdir
###        SP[CUBA.ROLLING_FRICTION] = '.true.' #tprnfor
        SP[CUBA.ROLLING_FRICTION] = 8 #ibrav
        SP[CUBA.ORIGINAL_POSITION] = [40,0.1,1.0] #write to celldm(1),(2),(3)
        #write n_atoms  nat - not needed since natoms can be derived from pc
        SP[CUBA.SCALING_COEFFICIENT] = 1 #1 atom type
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

        p1 = Particle([1.0* 1e-10,2.0* 1e-10,3.0* 1e-10])
        p1.data[CUBA.CHEMICAL_SPECIE] = 'C'
        p2 = Particle([2.0* 1e-10,3.0* 1e-10,4.0* 1e-10])
        p2.data[CUBA.CHEMICAL_SPECIE] = 'C'
        p3 = Particle([3.0* 1e-10,4.0* 1e-10,5.0* 1e-10])
        p3.data[CUBA.CHEMICAL_SPECIE] = 'C'
        p4 = Particle([4.0* 1e-10,5.0* 1e-10,6.0* 1e-10])
        p4.data[CUBA.CHEMICAL_SPECIE] = 'C'
        pc.add_particles([p1,p2,p3,p4])
        wrp.WriteEspressoInputFile(espresso_input_filename)

    def test_espresso_ppfile_write(self,ppfilename="testpp.in"):
        wrp = espresso_class.qe_functions()
        wrp.WriteEspressoPPFile(self,ppfilename="testpp.in"):



if __name__ == '__main__':
    unittest.main()