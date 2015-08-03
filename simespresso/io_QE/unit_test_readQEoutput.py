__author__ = 'jeremy'

import unittest
import espresso_input_file_read
import subprocess
import logging
import sys
from simphony.core.data_container import DataContainer
from simphony.core.cuba import CUBA
from simphony.core.keywords import KEYWORDS
from simphony.cuds.particles import Particle, Particles
import espresso_input_file_write

class OutcomesTest(unittest.TestCase):

    def setUp(self):
        self.filename = 'pwtest.in'
        with open(self.filename,'w') as text_file:
            text_file.write(_data_file_contents)

    def test_espresso_data_file_read(self):
        print('starting test of data file handler')
    #   What are the allowed keywords?
    #    for kw in KEYWORDS:
    #        print('kw:'+kw+'='+str(KEYWORDS[kw]))
        particle_container = espresso_input_file_read.ReadEspressoInputFile(self.filename)
        expected_atom_positions = []
        expected_atom_species = []
        expected_atom_positions.append((1.0e-10,2.0e-10,3.0e-10))
        expected_atom_positions.append((2.0e-10,3.0e-10,4.0e-10))
        expected_atom_positions.append((3.0e-10,4.0e-10,5.0e-10))
        expected_atom_positions.append((4.0e-10,5.0e-10,6.0e-10))
        expected_atom_species.append('C')
        expected_atom_species.append('C')
        expected_atom_species.append('C')
        expected_atom_species.append('C')
        i = 0

        #look for the above atoms in the particle container.
        #there is probably some intelligent way to do this but since the particle uid is
        #generated automatically (always, afaict) I can't find a way to index the particles as they appear in the original file.

        for particle in particle_container.iter_particles():
#            print('expected:'+str(expected_atom_positions[i])+' actual:'+str(particle))
#            print('data:'+str(particle.data))
            for expected_atom in expected_atom_positions:
 #               print('cur atom '+str(expected_atom)+' trying to match '+str(particle.coordinates))
                match = True
                for dimension in range(0,3):
                    if(expected_atom[dimension]!=particle.coordinates[dimension]):
                        match = False
 #                       print('no match:'+str(expected_atom[dimension])+'!='+str(particle.coordinates[dimension]))
                        break

                if match == True:
                    expected_atom_positions.remove(expected_atom)
                    self.assertTrue(particle.data[CUBA.CHEMICAL_SPECIE] == expected_atom_species[i])
  #                  print('cur atom matches')
                    break
        #we expect an empty list by now since every atom should have been matched and removed from list
        #of expect atoms
        self.assertTrue(expected_atom_positions==[])

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

