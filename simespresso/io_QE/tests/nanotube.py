__author__ = 'jeremy'

import os
import unittest

from simespresso.io_QE import espresso_class
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.particles import Particle, Particles


class nanotubes(object):

    def __init__(self):
        self.qe_wrapper = espresso_class.qe_functions()

    def place_atoms(self,length,diameter,helicity):
        self.qe_wrapper.pc = Particles('quantum_espresso_particles')
        particles = []
        for i in range(0,10):
            p = Particle([1.0 * 1e-10, 2.0 * 1e-10, 3.0 * 1e-10 * i])
            p.data[CUBA.CHEMICAL_SPECIE] = 'C'
            particles.append(p)
        self.qe_wrapper.pc.add_particles(particles)


        wrp.write_espresso_input_file(espresso_input_filename)



    def test_espresso_data_file_read(self):
        expected_atom_positions.append((1.0e-10, 2.0e-10, 3.0e-10))
        for particle in wrapper.pc.iter_particles():
                            particle.coordinates[dimension]):
                    self.assertTrue(particle.data[CUBA.CHEMICAL_SPECIE]


        wrp.pc = Particles('quantum_espresso_particles')




        wrapper.start_qe(name_in, name_out, path_to_espresso=path_to_espresso,
                         mpi=mpi, mpi_Nprocessors=mpi_Nprocessors)

if __name__ == '__main__':
    nt = nanotubes()
