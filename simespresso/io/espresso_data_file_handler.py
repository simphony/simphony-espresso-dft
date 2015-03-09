import string
from enum import Enum


class EspressoDataFileParser(object):
    """  This class parses  Espresso data files, either input or output
    (produced by the espresso command  write_data) and calls a handler 
    which processes the parsed information.
    A handler class is given the parsed information. This handler class can
    then determine what to  do with it.  For, example it could just store
    the data in memory (see LammpsSimpleDataHandler) or write it some other
    data file (e.g. a CUDS-file).

    cheatsheet for quantum espresso 
 
    0. create/obtain input file such as pp.in from cuds data
    1. run  mpi for instance 
            mpirun -np 48 /usr/local/espresso/bin/pp.x < input_pp.in > pp.out &
    2. convert output (which is a charge density file) into simphony format (see charge_density_xyz.cpp)
    

    Parameters
    ----------
    handler :
       handler will handle the parsed information provided by this class """
    def __init__(self, handler):
        self._handler = handler

    def parse(self, file_name):
        """ Read in data file containing starting state of simulation - 
            atom ids and positions, simulation type
        """
        self._handler.begin()
        state = _ReadState.UNKNOWN

        with open(file_name, 'r') as f:
            try:
                line_number = 0
                for line in f:
                    line_number += 1
                    #DEBUG
                    print('got line:'+str(line))
                    # skip blank lines
                    if not line.strip():
                        continue

                    state = _ReadState.get_state(state, line)
                    if state is _ReadState.ATOM_TYPES:
                        #DEBUG
                        print('reading atom types')
                        number_types = int(string.split(line, " ", 1)[0])
                        self._handler.process_number_atom_types(
                            number_types)
                    elif state is _ReadState.MASSES:
                        #DEBUG
                        print('reading masses')
                        values = line.split()
                        self._handler.process_masses(
                            int(values[0]),
                            int(values[1]))
                    elif state is _ReadState.ATOMS:
                        #DEBUG
                        print('reading atom ')
                        values = line.split()
                        id = int(values[0])
                        type_coord_etc = [int(values[1])]
                        for v in map(float, values[2:]):
                            type_coord_etc.append(v)
                        self._handler.process_atoms(id, type_coord_etc)
                    elif state is _ReadState.VELOCITIES:
                        values = line.split()
                        self._handler.process_velocities(
                            int(values[0]),
                            map(float, values[1:]))
                    else:
                        continue
            except Exception:
                print("problem with line number=", line_number, line)
                raise
        self._handler.end()


class _ReadState(Enum):
    UNKNOWN, UNSUPPORTED, \
        MASSES_BEGIN, MASSES, \
        VELOCITIES_BEGIN, VELOCITIES, \
        ATOM_TYPES, \
        ATOMS, \
        ATOMS_BEGIN = range(9)

    @staticmethod
    def get_state(current_state, line):
        """ Reads line and returns state
        """

        new_state = current_state

        # TODO how the state is determined and how
        # we transition to other states needs to be
        # rewritten.
        if (current_state is _ReadState.ATOMS_BEGIN or
                current_state is _ReadState.ATOMS):
            new_state = _ReadState.ATOMS
        elif (current_state is _ReadState.MASSES_BEGIN or
                current_state is _ReadState.MASSES):
            new_state = _ReadState.MASSES
        elif (current_state is _ReadState.VELOCITIES_BEGIN or
                current_state is _ReadState.VELOCITIES):
            new_state = _ReadState.VELOCITIES
        else:
            new_state = _ReadState.UNKNOWN

        if "atom types" in line:
            new_state = _ReadState.ATOM_TYPES
        elif "Masses" in line:
            new_state = _ReadState.MASSES_BEGIN
        elif "Velocities" in line:
            new_state = _ReadState.VELOCITIES_BEGIN
        elif "Atoms" in line:
            new_state = _ReadState.ATOMS_BEGIN
        elif "Pair Coeffs" in line:
            new_state = _ReadState.UNSUPPORTED
        elif new_state is _ReadState.UNKNOWN:
            new_state = _ReadState.UNSUPPORTED

        return new_state

