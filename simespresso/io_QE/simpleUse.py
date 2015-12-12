import numpy
import math
from simphony.core.cuba import CUBA
from simphony.core.cuds_item import CUDSItem
from simphony.cuds.particles import Particle, Particles
from simphony.io.h5_cuds import H5CUDS

#    from simphony.engine import quantumESPRESSO
# todo  figure out how to get my wrapper into simphony.engine
# instead of the following
from simespresso.io_QE import espresso_class

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
pc = Particles("Copper")
for base_vector in basis:
    position = [component*a_latt for component in base_vector]
    p = Particle(coordinates=position)
    p.data[CUBA.CHEMICAL_SPECIE] = ['Cu']  # this should be later an enum...
                                           # like, CUBA.CHEMICAL.ELEMENTS.Cu
    p.data[CUBA.MASS] = 63.546             # this is the atomic mass
    pc.add_particles([p])
    print position

super_cell = [[x*a_latt for x in v] for v in enumerate(unit_cell)]
pc.data_extension = {espresso_class.CUBAExtension.BOX_VECTORS: super_cell}

# later (D1.6) the BC should be part of the cuds: mycuds.BOX_VECTORS
# define the wrapper to use.

wrapper = espresso_class.qe_wrapper()
# Define the BC component of the SimPhoNy application model:
wrapper.BC_extension[quantumESPRESSO.CUBAExtension.BOX_FACES] = ["periodic",
                                                        "periodic",
                                                        "periodic"]

pc_w = wrapper.add_dataset(pc)

wrapper.SP_extension[quantumESPRESSO.CUBAExtension.PSEUDO_POTENTIAL] = 'Cu.pz-d-hgh.UPF'

# good for now, this is a standard pseudopotential, later we shall have a better way (actually we have it now but not implemented yet)

wrapper.CM_extension[quantumESPRESSO.CUBAExtension.K_POINT_SAMPLING_METHOD] = CUBA.Monkhorst-Pack

wrapper.CM_extension[quantumESPRESSO.CUBAExtension.K_POINT_SAMPLING] = [3, 3, 3 ]



wrapper.run()

# the wrapper would add the CUBA.TOTAL_ENERGY to the data of the pc within the value of the total energy from the output (log file) of QE

pc = wrapper.get_dataset("Copper")

pc.data_extension[quantumESPRESSO.CUBAExtension.TOTAL_ENERGY] # should print the tot eng