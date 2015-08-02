import string
from enum import Enum
from simphony.cuds.abstractlattice import ABCLattice
from simphony.core.data_container import DataContainer
from simphony.io.data_container_description import Record

from simphony.core.data_container import DataContainer
from simphony.core.cuba import CUBA



class EspressoInputDataFileParser(object):
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
       handler will handle the parsed information provided by this class """
    def __init__(self, handler):
        self._handler = handler
        self.atomtypes = ["C","H","He","N","O","Na","Mg"]

    def parse(self, file_name):
        """ Read in data file containing starting state of simulation - 
            atom ids and positions, simulation type
        """
#        self._handler.begin()
        self.dc = DataContainer()
        state = _ReadState.UNKNOWN

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
                    print('in main loop, got line:'+str(line))
                    # skip blank lines
     #               if not line.strip():
     #                   continue

                    state = _ReadState.get_state(state,line)
                    if state is _ReadState.CONTROL:
                        print('reading control section')
                        line = self.process_control(file_iter)
                        continue
                    elif state is _ReadState.SYSTEM:
                        print('reading system')
                        line = self.process_system(file_iter)
                        continue
                    elif state is _ReadState.ELECTRONS:
                        print('reading atoms')
                        line = self.process_electrons(file_iter)
                        continue
                    elif state is _ReadState.ATOMIC_SPECIES:
                        print('reading atomic species')
                        line = self.process_atomic_species(file_iter)
                        continue
                    elif state is _ReadState.K_POINTS:
                        print('reading k points')
                        values = line.split()
                        line = self.process_k_points(file_iter,mode=values[1])
                        continue
                    elif state is _ReadState.ATOMIC_POSITIONS:
                        print('reading atomic positions')
                        values = line.split()
                        line = self.process_atomic_positions(file_iter,units=values[1])
                        if line == 'EOF':
                            return
                        else:
                            continue

                    line = file_iter.next()
            except StopIteration:
                print('eof reached')
                return
            except Exception:
                print("problem with line number=", line_number, line)
                raise
#        self._handler.end()

    def process_control(self,f):
        print('processing control section')
        line = f.next()
        while _ReadState.get_state(_ReadState.CONTROL,line) == _ReadState.CONTROL:
            values = [x.strip() for x in line.split('=')]
            print('line in control section:'+str(line))
            if "calculation" in line:
                values = line.split('=')
                calculation = values[1]
            elif "restart_mode" in line:
                restart_mode = values[1]
            elif "pseudo_dir" in line:
                pseudo_dir = values[1]
            elif "prefix" in line:
                prefix = values[1]
            elif "tprnfor" in line:
                tprnfor = values[1]
            elif "max_seconds" in line:
                max_seconds = float(values[1])
            elif "outdir" in line:
                outdir = values[1]

            line = f.next()
        return line

    def process_system(self,f):
        print('processing system section')
        line = f.next()
        celldm=[0,0,0]
        while _ReadState.get_state(_ReadState.SYSTEM,line) == _ReadState.SYSTEM:
            values = [x.strip() for x in line.split('=')]
            print('line in control section:'+str(line))
            if "ibrav" in line:
                ibrav = int(values[1])
            elif "celldm(1)" in line:
                celldm[0] = float(values[1])
            elif "celldm(2)" in line:
                celldm[1] = float(values[1])
            elif "celldm(3)" in line:
                celldm[2] = float(values[1])
            elif "nat" in line:
                n_atoms = int(values[1])
            elif "ntyp" in line:
                n_atom_types = int(values[1])
            elif "ecutwfc" in line:
                ecutwfc = float(values[1])
            elif "ecutrho" in line:
                ecutrho = float(values[1]) #maybe int
            elif "input_dft" in line:
                input_dft = values[1]
            line = f.next()
        return line


    def process_electrons(self,f):
        print('processing eletrons section')
        line = f.next()
        while _ReadState.get_state(_ReadState.ELECTRONS,line) == _ReadState.ELECTRONS:
            values = [x.strip() for x in line.split('=')]
            print('line in electrons section:'+str(line))
            if "mixing_mode" in line:
                mixing_mode = values[1]
            elif "mixing_beta" in line:
                mixing_beta = float(values[1])
            elif "conv_thr" in line:
                conv_thr = values[1] #numbers like 1.0d-7 might have to be converted to float
            line = f.next()
        return line

    def process_atomic_species(self,f):
        print('processing atomic species section')
        line = f.next()
        while _ReadState.get_state(_ReadState.ATOMIC_SPECIES,line) == _ReadState.ATOMIC_SPECIES:
            values = line.split()
            print('line in atomic species section:'+str(line))
#            print('atomtypes:'+str(self.atomtypes))

            if len(values)>0:
                if values[0] in self.atomtypes:
                    print("atom type:"+values[0])
                    #self.dc(CHEMICAL_SPECIE = values[0])
                    DataContainer(CHEMICAL_SPECIE=values[0])
                    mass = float(values[1])
                    potential_file = values[2]
            line = f.next()
        return line

    def process_k_points(self,f,mode='automatic'):
        #skip line
        print('processing k_points section')
        line = f.next()
        while _ReadState.get_state(_ReadState.ATOMIC_SPECIES,line) == _ReadState.ATOMIC_SPECIES:
            values = [x.strip() for x in line.split('=')]
            print('line in k points section:'+str(line))
            K_points = values
            line = f.next()
        return line

    def process_atomic_positions(self,f,units='(angstrom)'):
        print('processing atomic_positions section')
        atom_positions = []
        try:
            line = f.next()
        except StopIteration:
            return('EOF')

        while _ReadState.get_state(_ReadState.ATOMIC_SPECIES,line) == _ReadState.ATOMIC_SPECIES:
            print('line in atomic positions section:'+str(line))
            values = [x.strip() for x in line.split('=')]
            if values[0] in self.atomtypes:
                atomtype = values[0]
                atom_pos[0] = values[1]
                atom_pos[1] = values[2]
                atom_pos[2] = values[3]
                atom_positions.append(atom_pos)
            try:
                line = f.next()
            except StopIteration:
                return('EOF')
        return line


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

        # TODO how the state is determined and how
        # we transition to other states needs to be
        # rewritten.
       # print('entering getstate')

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
        print('current state:'+str(new_state))
        return new_state

if __name__ == "__main__":
    filename = 'pw.in'
    print('started parsing file '+str(filename))
    parser = EspressoInputDataFileParser(filename)