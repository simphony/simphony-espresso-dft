from simphony.core.cuba import CUBA
from simphony.cuds.particles import Particle, Particles
from simespresso import qe_wrapper
import generate_nanotube


def generate_and_simulate_cswnt():
    particle_coordinates = generate_nanotube.read_xyz('tests/cswnt_m1n2l5.xyz')
    # particle_coordinates = generate_nanotube.generate_nanotube_xyz('armchair')

    unit_cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    a_latt = 3.61  # Angstroms
    super_cell = [[x*a_latt for x in v] for v in unit_cell]

    pc = Particles("Carbon")
    for particle in particle_coordinates:
        p = Particle(coordinates=particle)
        p.data[CUBA.CHEMICAL_SPECIE] = ['C']
        p.data[CUBA.MASS] = 12.0
        pc.add_particles([p])

    pc.data_extension = {'BOX_VECTORS': super_cell}
    wrapper = qe_wrapper.QeWrapper()
    wrapper.BC_extension['BOX_FACES'] = ["periodic", "periodic", "periodic"]
    wrapper.add_dataset(pc)
    wrapper.SP_extension['PSEUDO_POTENTIAL'] = '06-C.GGA.fhi.UPF'
    wrapper.CM_extension['K_POINT_SAMPLING_METHOD'] = "automatic"
    wrapper.CM_extension['K_POINT_SAMPLING'] = [3, 3, 3, 0, 0, 0]
    wrapper.CM_extension['DESIRED_SIMULATIONS'] = \
        ['TOTAL_ENERGY', 'CHARGE_DENSITY']
    wrapper.run()
    extracted_pc = wrapper.get_dataset("Carbon")
    charge_density = wrapper._data_manager._read_espresso_output_file('output.charge')
    print('checking particles:')
    for particle in extracted_pc.iter_particles():
        print('particle:'+str(particle))
    if 'TOTAL_ENERGY' in extracted_pc.data_extension:
        etot = extracted_pc.data_extension['TOTAL_ENERGY']  # should print eng
        print('tot energy:'+str(etot))
    else:
        print('tot energy not found')

if __name__ == "__main__":
    generate_and_simulate_cswnt()
