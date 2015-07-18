import string
from enum import Enum
from simphony.cuds.abstractlattice import ABCLattice

from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
#from simespresso.cuba_extension import CUBAExtension
from simphony.cuds.particles import Particle, Particles

#from simphony.cuds.abstractparticles import ABCParticles

import uuid


def ReadEspressoInputFile(file_name):
    """  This class parses  Espresso data files, either input or output
    (produced by the espresso command  write_data) and calls a handler 
    which processes the parsed information.
    A handler class is given the parsed information. This handler class can
    then determine what to  do with it.  For, example it could just store
    the data in memory (see LammpsSimpleDataHandler) or write it some other
    data file (e.g. a CUDS-file).

    cheatsheet for quantum espresso 
 
    0. create/obtain input file such as pp.in from cuds data
    1. run  (using mpi for instance )
        for file describing simulation
            mpirun -np 48 /usr/local/espresso/bin/pw.x < input_pw.in > pw.out &
        once simulation is done, run on file describing desired output
            mpirun -np 48 /usr/local/espresso/bin/pp.x < input_pp.in > pp.out &
    2. convert output (which is a charge density file) into simphony format (see charge_density_xyz.cpp)
    

    Parameters
    ----------
    handler :
       handler will handle the parsed information provided by this class


                    SP[CUBA.TORQUE] = calculation_type
            SP[CUBA.ZETA_POTENTIAL] = restart_mode
            SP[CUBA.YOUNG_MODULUS] = pseudo_dir
            SP[CUBA.VOLUME_FRACTION] = prefix
            SP[CUBA.AMPHIPHILICITY] = tprnfor
            SP[CUBA.NUMBER_OF_TIME_STEPS] = max_seconds
            SP[CUBA.OUTDIR] = outdir



        """


#CM.[CUBA.]
#BC.[CUBA.]

    state = _ReadState.UNKNOWN
#    BC = DataContainer()   #boundary conditions
#    CM = DataContainer()   #computational method
    SP = DataContainer()   #System Parameters and Conditions
#    SD = DataContainer()  #state data
    pc = Particles('quantum_espresso_particles')
    dc = DataContainer()
    with open(file_name, 'r') as f:
        line_number = 0
        file_iter = iter(f)
        line = file_iter.next()
        try:
            while line is not None:
                #file_iter.hasnext():
#                   #line = f[line_number]
                line_number += 1
                #DEBUG
                print('read line:'+str(line))
                # skip blank lines
 #               if not line.strip():
 #                   continue

                state = _ReadState.get_state(state,line)
                if state is _ReadState.CONTROL:
                    print('reading control section')
                    line = process_control(file_iter,SP)
                    continue
                elif state is _ReadState.SYSTEM:
                    print('reading system')
                    line = process_system(file_iter,SP)
                    continue
                elif state is _ReadState.ELECTRONS:
                    print('reading electrons')
                    line = process_electrons(file_iter)
                    continue
                elif state is _ReadState.ATOMIC_SPECIES:
                    print('reading atomic species')
                    line = process_atomic_species(file_iter)
                    continue
                elif state is _ReadState.K_POINTS:
                    print('reading k points')
                    values = line.split()
                    line = process_k_points(file_iter,mode=values[1])
                    continue
                elif state is _ReadState.ATOMIC_POSITIONS:
                    print('reading atomic positions')

                    values = line.split()
                    pc = process_atomic_positions(file_iter,pc,units=values[1])
                    break

                line = file_iter.next()
        except StopIteration:
            print('eof reached')
        except Exception:
            print("problem with line number=", line_number, line)
            raise
    return pc


def process_control(f,SP):
    print('processing control section')
    line = f.next()
    while _ReadState.get_state(_ReadState.CONTROL,line) == _ReadState.CONTROL:
        values = [x.strip() for x in line.split('=')]
        print('line in control section:'+str(line))
        if "calculation" in line:
            values = line.split('=')
            calculation_type = values[1]
            #THIS IS A HACK . Use of ZETA POTENTIAL for calculation_type is because attempts to use
            #CUBAExtension were unsuccessful.
            SP[CUBA.TORQUE] = calculation_type
        elif "restart_mode" in line:
            restart_mode = values[1]
            #THIS IS A HACK . Use of ZETA POTENTIAL for restart mode
            SP[CUBA.ZETA_POTENTIAL] = restart_mode
        elif "pseudo_dir" in line:
            pseudo_dir = values[1]
            #THIS IS A HACK . Use of YOUNG MODULUS for pseudo_dir
            SP[CUBA.YOUNG_MODULUS] = pseudo_dir
        elif "prefix" in line:
            prefix = values[1]
            #THIS IS A HACK . Use of VOLUME FRACTION for prefix
            SP[CUBA.VOLUME_FRACTION] = prefix
        elif "tprnfor" in line:  #calculate forces
            tprnfor = values[1]
            #THIS IS A HACK . Using AMPHILICITY for tprnfor
            SP[CUBA.AMPHIPHILICITY] = tprnfor
        elif "max_seconds" in line:
            max_seconds = float(values[1])
            #THIS IS A HACK . Using NUMBER_OF_TIME_STEPS for max_seconds
            SP[CUBA.NUMBER_OF_TIME_STEPS] = max_seconds
        elif "outdir" in line:
            outdir = values[1]
            #THIS IS A HACK . Using DIRECTION for outdir
            SP[CUBA.DIRECTION] = outdir

        line = f.next()
    return line

def process_system(f,SP):
    print('processing system section')
    line = f.next()
    celldm=[0,0,0]
    while _ReadState.get_state(_ReadState.SYSTEM,line) == _ReadState.SYSTEM:
        values = [x.strip() for x in line.split('=')]
        print('line in control section:'+str(line))
        if "ibrav" in line:  #bravais lattice index
            ibrav = int(values[1])
#            SP[CUBAExtension.IBRAV] = ibrav
        elif "celldm(1)" in line:
            celldm[0] = float(values[1])
 #           SP[CUBA.LATTICE_VECTORS = celldm[0]
        elif "celldm(2)" in line:
            celldm[1] = float(values[1])
#            SP[CUBA.BOX_VECTORS][1] = celldm[1]
        elif "celldm(3)" in line:
            celldm[2] = float(values[1])
#            SP[CUBA.BOX_VECTORS][2] = celldm[2]
            SP[CUBA.LATTICE_VECTORS] = celldm

        elif "nat" in line:
            n_atoms = int(values[1])
        elif "ntyp" in line:
            n_atom_types = int(values[1])
        elif "ecutwfc" in line:
            ecutwfc = float(values[1])
#            SP[CUBAExtension.KINETIC_ENERGY_CUTOFF_FOR_WAVEFUNCTIONS] = ecutwfc
        elif "ecutrho" in line:
            ecutrho = float(values[1]) #maybe int
#            SP[CUBAExtension.KINETIC_ENERGY_CUTOFF_FOR_CHARGE_DENSITY_AND_POTENTIAL] = ecutwfc
        elif "input_dft" in line:
            input_dft = values[1]
#            SP[CUBAExtension.EXCHANGE_CORRELATION_FUNCTIONAL] = input_dft
        line = f.next()
    return line


def process_electrons(f):
    print('processing eletrons section')
    line = f.next()
    while _ReadState.get_state(_ReadState.ELECTRONS,line) == _ReadState.ELECTRONS:
        values = [x.strip() for x in line.split('=')]
        print('line in electrons section:'+str(line))
        if "mixing_mode" in line:
            mixing_mode = values[1]
#            SP[CUBAExtension.MIXING_MODE] = mixing_mode
        elif "mixing_beta" in line:
            mixing_beta = float(values[1])
#            SP[CUBAExtension.MIXING_BETA] = mixing_beta
        elif "conv_thr" in line:
            conv_thr = values[1] #numbers like 1.0d-7 might have to be converted to float
#            SP[CUBAExtension.CONVERGENCE_THRESHOLD] = conv_thr
        line = f.next()
    return line

def process_atomic_species(f):
    print('processing atomic species section')
    line = f.next()
    while _ReadState.get_state(_ReadState.ATOMIC_SPECIES,line) == _ReadState.ATOMIC_SPECIES:
        values = line.split()
        print('line in atomic species section:'+str(line))
#            print('atomtypes:'+str(self.atomtypes))

        if len(values)>0:
            if values[0] in atomtypes:
                print("atom type:"+values[0])
                #self.dc(CHEMICAL_SPECIE = values[0])
#                DataContainer(CHEMICAL_SPECIE=values[0])
                mass = float(values[1])
                potential_file = values[2]
        line = f.next()
    return line

def process_k_points(f,mode='automatic'):
    #skip line
    print('processing k_points section')
    line = f.next()
#    SP[CUBAExtension.K_POINTS_MODE] = mode
    i = 0
    while _ReadState.get_state(_ReadState.K_POINTS,line) == _ReadState.K_POINTS:
#        print('line:'+str(line))
        values = line.split()
        if len(values):
            K_points = (values)
            print('k points:'+str(K_points))
#            SP[CUBAExtension.K_POINTS] = K_points

        line = f.next()
    return line


#I am not sure whether to use datacontainer or particlecontainer - maybe PC goes in DC?
def process_atomic_positions(f,pc,units='(angstrom)'):
    print('processing atomic_positions section')
    try:
        line = f.next()
    except StopIteration:
        return('EOF')

    while _ReadState.get_state(_ReadState.ATOMIC_POSITIONS,line) == _ReadState.ATOMIC_POSITIONS:
        print('line in atomic positions section:'+str(line))
        values = line.split()
        print('values:'+str(values))
        atom_pos = [0,0,0]
        i = 0
        if values[0] in atomtypes:
            atomtype = values[0]
            atom_id = i
            atom_pos[0] = float(values[1]) * 1e-10  #store position in meters; original in Angstrom
            atom_pos[1] = float(values[2]) * 1e-10
            atom_pos[2] = float(values[3]) * 1e-10
            s = str(i)
            #make uid using string of index
            #uidstring =
            #uid = uuid.UUID(uidstring)

            p = Particle([atom_pos[0],atom_pos[1],atom_pos[2]])  #,uuid.UUID(int=i)
            #from cuds/tests/test_particles.py: particle = Particle([20.5, 30.5, 40.5], uuid.UUID(int=33), data)
            print('uid:'+str(p.uid))
            p.data[CUBA.CHEMICAL_SPECIE] = atomtype

            part = Particle()
            part_container = Particles(name="foo")
            part_container.add_particle(part)


            pc.add_particle(p)
  #          print('pc:'+str(pc))
            i = i +1
        try:
            line = f.next()
        except StopIteration:
            return('EOF')

    return pc


class _ReadState(Enum):
    UNKNOWN, UNSUPPORTED, \
        CONTROL,\
        SYSTEM, \
        ELECTRONS, \
        IONS, \
        CELL, \
        ATOMIC_SPECIES, \
        K_POINTS, \
        ATOMIC_POSITIONS = range(10)

    @staticmethod
    def get_state(current_state, line):
        """ Reads line and returns state
        """
        new_state = current_state

        if "&CONTROL" in line:
            new_state = _ReadState.CONTROL
        elif "&SYSTEM" in line:
            new_state = _ReadState.SYSTEM
        elif "&ELECTRONS" in line:
            new_state = _ReadState.ELECTRONS
        elif "&IONS" in line:
            new_state = _ReadState.UNSUPPORTED
        elif "&CELL" in line:
            new_state = _ReadState.UNSUPPORTED
        elif "ATOMIC_SPECIES" in line:
            new_state = _ReadState.ATOMIC_SPECIES
        elif "K_POINTS" in line:
            new_state = _ReadState.K_POINTS
        elif "ATOMIC_POSITIONS" in line:
            new_state = _ReadState.ATOMIC_POSITIONS
  #      print('current state:'+str(new_state))
        return new_state


atomtypes = ["C","H","He","N","O","Na","Mg"]

if __name__ == "__main__":
    filename = '../../examples/input_pw.in'
    print('started parsing file '+str(filename))
    ReadEspressoInputFile(filename)
