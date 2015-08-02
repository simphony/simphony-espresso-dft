import string
from enum import Enum
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.particles import Particle, Particles


def WriteEspressoInputFile(file_name,data_container):
    """
    :param file_name: name of the input file to write
    :param data_container: contains data to write to file
    :return:
    """
    d=DataContainer.items()
    pc=data_container
    for item in d:

        with open(file_name, 'w') as f:
            try:

                f.write()
            except:
                ('error in write block of WriteEspressoInputFile')
                raise
    f.closed



def process_control(f):
    print('processing control section')
    line = f.next()
    while _ReadState.get_state(_ReadState.CONTROL,line) == _ReadState.CONTROL:
        values = [x.strip() for x in line.split('=')]
        print('line in control section:'+str(line))
        if "calculation" in line:
            values = line.split('=')
            calculation_type = values[1]
            #THIS IS A HACK . Use of ZETA POTENTIAL for restart mode is because attempts to use
            #CUBAExtension were unsuccessful.
            SP[CUBA.TORQUE] = calculation_type
        elif "restart_mode" in line:
            restart_mode = values[1]
            #THIS IS A HACK . Use of ZETA POTENTIAL for restart mode is because attempts to use
            #CUBAExtension were unsuccessful.
            SP[CUBA.ZETA_POTENTIAL] = restart_mode
        elif "pseudo_dir" in line:
            pseudo_dir = values[1]
            #THIS IS A HACK . Use of YOUNG MODULUS for restart mode is because attempts to use
            #CUBAExtension were unsuccessful.
            SP[CUBA.YOUNG_MODULUS] = pseudo_dir
        elif "prefix" in line:
            prefix = values[1]
            #THIS IS A HACK . Use of VOLUME FRACTION for prefix is because attempts to use
            #CUBAExtension were unsuccessful.
            SP[CUBA.VOLUME_FRACTION] = prefix
        elif "tprnfor" in line:  #calculate forces
            tprnfor = values[1]
#            SP[CUBAExtension.TPRNFOR] = tprnfor
        elif "max_seconds" in line:
            max_seconds = float(values[1])
#            SP[CUBAExtension.SIMULATION_MAXIMUM_RUNTIME] = max_seconds
        elif "outdir" in line:
            outdir = values[1]
#            SP[CUBAExtension.OUTDIR] = outdir

        line = f.next()
    return line

def process_system(f):
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
#            SP[CUBAExtension.BOX_VECTORS][0] = celldm[0]
        elif "celldm(2)" in line:
            celldm[1] = float(values[1])
#            SP[CUBAExtension.BOX_VECTORS][1] = celldm[1]
        elif "celldm(3)" in line:
            celldm[2] = float(values[1])
#            SP[CUBAExtension.BOX_VECTORS][2] = celldm[2]
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
            atom_pos[0] = float(values[1]) * 1e-10  #position in meters, original in Angstrom
            atom_pos[1] = float(values[2]) * 1e-10
            atom_pos[2] = float(values[3]) * 1e-10
            s = str(i)
            #make uid using string of index
            #uidstring =
            #uid = uuid.UUID(uidstring)

            p = Particle(coordinates=[atom_pos[0],atom_pos[1],atom_pos[2]])

            print('uid:'+str(p.uid))
            p.data[CUBA.CHEMICAL_SPECIE] = atomtype
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
