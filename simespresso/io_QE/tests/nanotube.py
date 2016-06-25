__author__ = 'jeremy'
from simespresso.io_QE import qe_file_io
from simphony.core.cuba import CUBA
from simphony.cuds.particles import Particle, Particles


class nanotubes(object):

    def __init__(self):
        self.qe_wrapper = qe_file_io.qe_functions()

    def make_and_test_nanotube(self):
        length = 10**-9
        diameter = 1**-10
        helicity = 2
        self.qe_wrapper.pc = Particles('quantum_espresso_particles')
        self.place_atoms(self, pc, length, diameter, helicity)

    def place_atoms(self,length,diameter,helicity):
        particles = []
        for i in range(0,10):
#vjust testing this is obviously not a nanotube
            p = Particle([1.0 * 1e-10, 2.0 * 1e-10, 3.0 * 1e-10 * i])
            p.data[CUBA.CHEMICAL_SPECIE] = 'C'
            particles.append(p)
        pc = self.qe_wrapper.pc
        pc.add_particles(particles)

    def write_qefiles(self):
        pass
        # wrp.write_espresso_input_file(espresso_input_filename)
        # particle.coordinates[dimension]):
        # particle.data[CUBA.CHEMICAL_SPECIE]
        # wrp.pc = Particles('quantum_espresso_particles')

        # wrapper.start_qe(name_in, name_out, path_to_espresso=path_to_espresso,
        #                 mpi=mpi, mpi_Nprocessors=mpi_Nprocessors)

if __name__ == '__main__':
    nt = nanotubes()
    nanotubes.make_and_test_nanotube()