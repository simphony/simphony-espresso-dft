

    import numpy

    import math

    from simphony.core.cuds import CUDS

    from simphony.core.cuba import CUBA

    from simphony.core.cuds_item import CUDSItem

    from simphony.cuds.particles import Particle, Particles

    from simphony.engine import quantumESPRESSO

    from simphony.io.h5_cuds import H5CUDS



################

### this part will be an utility later on...

# Create the Cu unit cell, assuming a simple cubic system with 4 basis atoms (for an FCC latticle).



# The lattice parameter (in a cubic setup)

a_latt = 3.61e  # this is in Angstroms, Please assume Angstroms, and atomic units for the rest, and atomic mass unit, later we shall make it uniform.



# Use a SC unit cell with basis for the FCC system

unit_cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

# The basis of the FCC system in the SC setup:

basis = [

    [0.0, 0.0, 0.0],

    [0.5, 0.5, 0.0],

    [0.5, 0.0, 0.5],

    [0.0, 0.5, 0.5]

]

###



# create an empty cuds, a common universal data structure to hole all the data relevant for a

# computational model, including solver parameters (SP), Physics Equations (PE) and Material Relations, as well as the

# data sets.



# define the top model Physics Equations:

pe = DFT()

pe.electronc = DFT.pseudopotential() # it does not take any parameters now, it only says that we use pseudopotentials, we shall extend this with time to have a type attribute specifying the differetnt flavours of PP.



# create an empty CUDS:

cuds = CUDS(name = "CuTotEng", description="calculate total energy of Cu from DFT")



# add the pe:

cuds.add_pe(pe)



# define the material and its common properties:

copper = Material(name = "copper", description="Cu unit cell for testing purposes", data=DataContainer(CUBAL.LATTICE_PARAMETER:alat))



# add the material to the cuds, there can be more than one

cuds.add_material(copper)



# define the dataset

pc = Particles("Copper")

cuds = CUDS ("Copper DFT Simulation")

for pos in basis:

    pos2 = [pos[0]*a_latt, pos[1]*a_latt, pos[2]*a_latt]

    p = Particle(coordinates=pos2)

    p.data[CUBA.CHEMICAL_SPECIE] = ['Cu']  # this should be later an enum...

                                           # like, CUBA.CHEMICAL.ELEMENTS.Cu

    p.data[CUBA.MASS] = 63.546             # this is the atomic mass

    p.data[CUBA.MATERIAL] = copper.uuid # this tells to which material each point

                                        # belongs too

    pc.add_particles([p])



CUDS.add_dataset(pc)

# the cuds collects all datasets and manages them internally. The most straightforward thing is to establish a mapping using uuids of the pc and internal uuid of the engine, and then save the atom positions in the QE file, at any point the cuds can then fetch it back as needed. No need to manage two buffers!



super_cell = [

    tuple(x*a_latt for x in v) for i, v in enumerate(unit_cell)]



cuds.box = super_cell # the cuds manages the supercell too.





# Define the BC component of the SimPhoNy application model:

bc = BoxBoundaryCondition(CUBA.PERRIODIC, CUBA.PERIODIC, CUBA.PERIODIC)

cuds.bc = BC # cuds also manages the Box type of boundary conditions



pp = UPF_PseudoPotential(name="PseudoPotentialCu", description="pseudopotential of Cu, it contains default cutoffs etc, so no need to specify them", materials=[copper.uuid], pseudopotential_name = 'Cu.pz-d-hgh.UPF' )



cuds.add_material_relation (pp)



kpoints = monkhorst_pack (grid=[3, 3, 3 ])

cuds.add_solver_par(kpoints)

# this would be managed int he SP (solver parameters, note the name change) of the cuds.





# now we could save the cuds to h5cuds for later but now we want to calculate since we are impatient...



# define the wrapper to use.

wrapper = quantumESPRESSO.quantumESPRESSOWrapper()



# now send the cuds computational model data structure to the wrapper:

wrapper.cuds = cuds

# or wrapper.set_cuds (cuds)



# run the wrapper:

out = wrapper.run()  # or wrapper.run(cuds)



# now we need to get the energy back!

print out.total_energy