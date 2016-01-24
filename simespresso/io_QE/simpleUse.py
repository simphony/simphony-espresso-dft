from simphony.core.cuba import CUBA
from simphony.cuds.particles import Particle, Particles
import os
import logging

#    from simphony.engine import quantumESPRESSO
# todo  figure out how to get wrapper into simphony.engine
# instead of the following

from simespresso import qe_wrapper
#import qeCubaExtensions
#import espresso_class


from qeCubaExtensions import qeCUBAExtension

# Create the Cu unit cell, assuming a simple cubic system with 4 basis
# atoms (for an FCC latticle)
# The lattice parameter (in a cubic setup)
a_latt = 3.61  # this is in Angstroms, Please assume Angstroms, and atomic
# units for the rest, and atomic mass unit, later we shall make it uniform.
# Use a SC unit cell with basis for the FCC system

unit_cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
# The basis of the FCC system in the SC setup:
basis = [
    [0.0, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5] ]
basis = [
    [0.0, 0.0, 0.0],
    [0.0, 0.5, 0.5]]
pc = Particles("Copper")
for base_vector in basis:
    position = [component * a_latt for component in base_vector]
    p = Particle(coordinates=position)
    p.data[CUBA.CHEMICAL_SPECIE] = ['Cu']  # this should be later an enum...
                                           # like, CUBA.CHEMICAL.ELEMENTS.Cu
    p.data[CUBA.MASS] = 63.546             # this is the atomic mass
    pc.add_particles([p])
    logging.debug('position:'+str(position))
super_cell = [[x*a_latt for x in v] for v in unit_cell]
pc.data_extension = {qeCUBAExtension.BOX_VECTORS: super_cell}

# later (D1.6) the BC should be part of the cuds: mycuds.BOX_VECTORS
# define the wrapper to use.

wrapper = qe_wrapper.QeWrapper()
# Define the BC component of the SimPhoNy application model:
#wrapper.BC_extension[qeCUBAExtension.BOX_FACES] = ["periodic",
 #                                                       "periodic",
  #                                                      "periodic"]
wrapper.BC_extension['BOX_FACES'] = ["periodic",
                                                        "periodic",
                                                        "periodic"]
pc_w = wrapper.add_dataset(pc)
#wrapper.SP_extension[qeCUBAExtension.PSEUDO_POTENTIAL] = 'Cu.pz-d-hgh.UPF'
cwd = os.path.dirname(__file__)
wrapper.SP_extension['PSEUDO_POTENTIAL'] = os.path.join(cwd,'Cu.pz-d-hgh.UPF')
wrapper.SP_extension['PSEUDO_POTENTIAL'] = os.path.join(cwd,'Cu.pz-d-rrkjus.UPF')
wrapper.SP_extension['PSEUDO_POTENTIAL'] = 'Cu.pz-d-rrkjus.UPF'

logging.debug('pp file:'+str(wrapper.SP_extension['PSEUDO_POTENTIAL']))
#internal = wrapper._data_manager.SP_extension['PSEUDO_POTENTIAL']
#logging.debug('pp file:'+str(internal))

# good for now, this is a standard pseudopotential,
# later we shall have a better way (actually we have it now
# but not implemented yet)
#wrapper.CM_extension['K_POINT_SAMPLING_METHOD'] = "Monkhorst-Pack"
wrapper.CM_extension['K_POINT_SAMPLING_METHOD'] = "automatic"
wrapper.CM_extension['K_POINT_SAMPLING'] = [3, 3, 3, 0, 0, 0]
wrapper.run()
# the wrapper would add the CUBA.TOTAL_ENERGY to the data of the pc within the value of the total energy from the output (log file) of QE
names = wrapper.get_dataset_names()
print('dataset names:'+str(names))
extracted_pc = wrapper.get_dataset("Copper")
print('checking particles:')
for particle in extracted_pc:
    print('particle:'+str(particle))
extracted_pc.data_extension[qeCUBAExtension.TOTAL_ENERGY] # should print the tot eng
