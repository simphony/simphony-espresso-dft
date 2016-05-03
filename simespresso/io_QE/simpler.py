from simphony.core.cuba import CUBA
from simphony.cuds.particles import Particle, Particles
import logging
import  matplotlib.pyplot as plt
from simespresso import qe_wrapper
from qeCubaExtensions import qeCUBAExtension

a_latt = 3.61  # this is in Angstroms, Please assume Angstroms, and atomic
unit_cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
basis = [
    [0.0, 0.0, 0.0],
    [0.0, 1.0, 0]]
pc = Particles("Copper")
for base_vector in basis:
    position = [component * a_latt for component in base_vector]
    p = Particle(coordinates=position)
    p.data[CUBA.CHEMICAL_SPECIE] = ['Cu']
    p.data[CUBA.MASS] = 63.546
    pc.add_particles([p])
    logging.debug('position:'+str(position))
super_cell = [[x*a_latt for x in v] for v in unit_cell]
pc.data_extension = {qeCUBAExtension.BOX_VECTORS: super_cell}
wrapper = qe_wrapper.QeWrapper()
wrapper.BC_extension['BOX_FACES'] = ["periodic",
                                                        "periodic",
                                                        "periodic"]
wrapper.add_dataset(pc)
wrapper.SP_extension['PSEUDO_POTENTIAL'] ='Cu.pz-d-hgh.UPF'
wrapper.CM_extension['K_POINT_SAMPLING_METHOD'] = "automatic"
wrapper.CM_extension['K_POINT_SAMPLING'] = [3, 3, 3, 0, 0, 0]
wrapper.run()
extracted_pc = wrapper.get_dataset("Copper")
print('checking particles:')
for particle in extracted_pc.iter_particles():
    print('particle:'+str(particle))
if 'TOTAL_ENERGY' in extracted_pc.data_extension:
    etot = extracted_pc.data_extension['TOTAL_ENERGY'] # should print the tot eng
    print('tot energy:'+str(etot))
else:
    print('tot energy not found')

