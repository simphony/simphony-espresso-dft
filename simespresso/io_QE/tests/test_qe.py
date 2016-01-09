__author__ = 'jeremy'

import os
import unittest

from simespresso.io_QE import qe_file_io
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.particles import Particle, Particles


class OutcomesTest(unittest.TestCase):

    def setUp(self):
        self.filename = 'pwtest.in'
        with open(self.filename, 'w') as text_file:
            text_file.write(_data_file_contents)

    def test_espresso_data_file_read(self):
        print('TEST OF READING QE INPUT FILE')
        wrapper = qe_file_io.qe_functions()
        wrapper.read_espresso_input_file(self.filename)
        expected_atom_positions = []
        expected_atom_species = []
        expected_atom_positions.append((1.0e-10, 2.0e-10, 3.0e-10))
        expected_atom_positions.append((2.0e-10, 3.0e-10, 4.0e-10))
        expected_atom_positions.append((3.0e-10, 4.0e-10, 5.0e-10))
        expected_atom_positions.append((4.0e-10, 5.0e-10, 6.0e-10))
        expected_atom_species.append('C')
        expected_atom_species.append('C')
        expected_atom_species.append('C')
        expected_atom_species.append('C')
        i = 0

        for particle in wrapper.pc.iter_particles():
            # print('expected:'+str(expected_atom_positions[i])+
            # ' actual:'+str(particle))
            # print('data:'+str(particle.data))
            for expected_atom in expected_atom_positions:
                # print('cur atom '+str(expected_atom)+
                # ' trying to match '+str(particle.coordinates))
                match = True
                for dimension in range(0, 3):
                    if(expected_atom[dimension] !=
                            particle.coordinates[dimension]):
                        match = False
                        # print('no match:'+str(expected_atom[dimension])+
                        # '!='+str(particle.coordinates[dimension]))
                        break

                if match is True:
                    expected_atom_positions.remove(expected_atom)
                    self.assertTrue(particle.data[CUBA.CHEMICAL_SPECIE]
                                    == expected_atom_species[i])
                    # print('cur atom matches')
                    break

    def test_espresso_data_file_write(self):
        print('TEST OF WRITING QE INPUT FILE')
        wrp = qe_file_io.qe_functions()
        espresso_input_filename = 'pw_generated.in'
        wrp.SP = DataContainer()
        wrp.pc = Particles('quantum_espresso_particles')
        # write parameters for a particular working input file
        SP = wrp.SP
        pc = wrp.pc
        SP[CUBA.TORQUE] = 'scf'
        # calculation
        SP[CUBA.ZETA_POTENTIAL] = 'from_scratch'
        # restart_mode
        SP[CUBA.YOUNG_MODULUS] = './'
        # pseudo dir
        SP[CUBA.VOLUME_FRACTION] = 'NANO'
        # prefix
        SP[CUBA.AMPHIPHILICITY] = '.true.'
        # tprnfor
        SP[CUBA.NUMBER_OF_TIME_STEPS] = 82000
        # max_seconds
        SP[CUBA.DIRECTION] = './'
        # SP[CUBA.ROLLING_FRICTION] = '.true.' #tprnfor
        SP[CUBA.ROLLING_FRICTION] = 8
        # ibrav
        SP[CUBA.ORIGINAL_POSITION] = [40, 0.1, 1.0]
        # write to celldm(1),(2),(3)
        # write n_atoms  nat - not needed since natoms can be derived from pc
        SP[CUBA.SCALING_COEFFICIENT] = 1
        # 1 atom type
        SP[CUBA.LN_OF_RESTITUTION_COEFFICIENT] = 60.0
        # ecutwfc
        SP[CUBA.POISSON_RATIO] = 240
        # ecutrho
        SP[CUBA.LATTICE_SPACING] = 'vdw-df-c09'
        # input_dft
        SP[CUBA.SMOOTHING_LENGTH] = 'local-TF'
        # mixing_mode
        SP[CUBA.PHASE_INTERACTION_STRENGTH] = 0.8
        # mixing_beta
        SP[CUBA.DEBYE_LENGTH] = '1.0d-7'
        # conv_thr
        SP[CUBA.CHEMICAL_SPECIE] = ['C']
        SP[CUBA.MASS] = [12.0107]
        SP[CUBA.FRICTION_COEFFICIENT] = ['06-C.GGA.fhi.UPF']
        SP[CUBA.PROBABILITY_COEFFICIENT] = 'automatic'  # mode
        SP[CUBA.EQUATION_OF_STATE_COEFFICIENT] = [1, 4, 1, 0, 0, 0]  # K_POINTS
        SP[CUBA.KINEMATIC_VISCOSITY] = '(angstrom)'  # ATOMIC_POSITIONS

        p1 = Particle([1.0 * 1e-10, 2.0 * 1e-10, 3.0 * 1e-10])
        p1.data[CUBA.CHEMICAL_SPECIE] = 'C'
        p2 = Particle([2.0 * 1e-10, 3.0 * 1e-10, 4.0 * 1e-10])
        p2.data[CUBA.CHEMICAL_SPECIE] = 'C'
        p3 = Particle([3.0 * 1e-10, 4.0 * 1e-10, 5.0 * 1e-10])
        p3.data[CUBA.CHEMICAL_SPECIE] = 'C'
        p4 = Particle([4.0 * 1e-10, 5.0 * 1e-10, 6.0 * 1e-10])
        p4.data[CUBA.CHEMICAL_SPECIE] = 'C'
        pc.add_particles([p1, p2, p3, p4])
        wrp.write_espresso_input_file(espresso_input_filename)

    def test_espresso_ppfile_write(self, ppfilename="testpp.in"):
        print('TEST OF WRITING QE PP FILE')

        wrp = qe_file_io.qe_functions()
        wrp.write_espresso_pp_file(ppfilename="testpp.in")

    def test_read_espresso_output_file(self):
        print('TEST OF READING QE OUTPUT FILE')

        file_name = 'pwtest.out'
        if not(os.path.exists(file_name)):
            # logging.debug("file "+str(file_name)+" not found")
            print("file "+str(file_name)+" not found")
            return(1)
        print('testing reading of qe output file '+str(file_name))
        qe_wrapper = qe_file_io.qe_functions()
        qe_wrapper.read_espresso_output_file(file_name)

    def test_running_index_to_node_inredex(self):
        print('testing espresso_data_file_read')
        n_latticepoints = [10, 7, 6]
        index = 5
        espresso_wrapper = qe_file_io.qe_functions()
        indices = espresso_wrapper.running_index_to_node_index(index,
                                                               n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices == [5, 0, 0])

        index = 15
        indices = espresso_wrapper.running_index_to_node_index(index,
                                                               n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices == [5, 1, 0])

        index = 23
        indices = espresso_wrapper.running_index_to_node_index(index,
                                                               n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices == [3, 2, 0])

        index = 70
        indices = espresso_wrapper.running_index_to_node_index(index,
                                                               n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices == [0, 0, 1])

        index = 75
        indices = espresso_wrapper.running_index_to_node_index(index,
                                                               n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices == [5, 0, 1])

        index = 150
        indices = espresso_wrapper.running_index_to_node_index(index,
                                                               n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices == [0, 1, 2])

    def test_start_qe(self):
        print('TEST OF STARTING QE')
        name_in = './pwtest.in'
        name_out = './pwtest.out'
        mpi = False
        mpi_Nprocessors = 2
        path_to_espresso = 'pw.x'  # /usr/bin/pw.x is not default apparently
        wrapper = qe_file_io.qe_functions()
#        wrapper.start_qe2(name_in, name_out)
        wrapper.start_qe(name_in, name_out, path_to_espresso=path_to_espresso,
                         mpi=mpi, mpi_Nprocessors=mpi_Nprocessors)

_data_file_contents = """&CONTROL
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = './',
    prefix='NANO'
    tprnfor = .true.
    max_seconds=82000
    outdir = './'
/
&SYSTEM
    ibrav=8
    celldm(1)=40
    celldm(2)= 0.115981
    celldm(3)=1.00
    nat=6
    ntyp=1
    ecutwfc=60.0
    ecutrho=240
    input_dft='vdw-df-c09'
 /
 &ELECTRONS
    mixing_mode='local-TF'
    mixing_beta = 0.8
    conv_thr =  1.0d-7
 /
 &IONS
 /
 &CELL
 /
ATOMIC_SPECIES
C 12.0107 06-C.GGA.fhi.UPF

K_POINTS automatic
  1 4 1 0 0 0

ATOMIC_POSITIONS (angstrom)
C        1.0  2.0   3.0
C        2.0  3.0   4.0
C        3.0  4.0   5.0
C        4.0  5.0   6.0
"""

if __name__ == '__main__':
    unittest.main()
