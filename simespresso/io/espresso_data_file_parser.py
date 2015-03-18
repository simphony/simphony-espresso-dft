import string
from enum import Enum


class EspressoDataFileParser(object):
    """  
    This class parses  Espresso input data files 
    (produced by the espresso command  write_data) and calls a handler 
    which processes the parsed information.
    A handler class is given the parsed information. This handler class can
    then store the data in memory, write it to a CUDS-file, etc.

    cheatsheet for quantum espresso 
 
    0. create/obtain input file such as pw.in from cuds data
    1. run  mpi for instance 
            mpirun -np 48 /usr/local/espresso/bin/pw.x < input_pw.in > pw.out &
    2. create output file of charge density  
            mpirun -np x /usr/local/espresso/bin/pp.x $ <$ name.in $ >$ name.out
    2. convert output (which is a charge density file) into simphony format (see charge_density_xyz.cpp)
    

    Parameters
    ----------
    handler :
       handler will handle the parsed information provided by this class 
    """
    
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
   	    	 	#DEBUG
		        print('current state:'+str(state))
                    if state is _ReadState.CONTROL:
			values = line.split()
			if "calculation" in line:
				calculation = values[1]			
			else if "restart_mode" in line:
				restart_mode = values[1]			
			else if "pseudo_dir" in line:
				pseudo_dir = values[1]			
			else if "prefix" in line:
				prefix = values[1]			
			else if "tprnfor" in line:
				tprnfor = values[1]			
			else if "max_seconds" in line:
				max_seconds = float(values[1])			
			else if "outdir" in line:
				outdir = values[1]			

                    elif state is _ReadState.SYSTEM:
                        values = line.split()
			if "ibrav" in line:
				ibrav = int(values[1])			
			else if "celldm(1)" in line:
				celldm[0] = float(values[1])			
			else if "celldm(2)" in line:
				celldm[1] = float(values[1])			
			else if "celldm(3)" in line:
				celldm[2] = float(values[1])			
			else if "nat" in line:
				n_atoms = int(values[1])			
			else if "ntyp" in line:
				n_atom_types = int(values[1])			
			else if "ecutwfc" in line:
				ecutwfc = float(values[1])			
			else if "ecutrho" in line:
				ecutrho = float(values[1]) //maybe int			
			else if "input_dft" in line:
				input_dft = values[1]			

                    elif state is _ReadState.ELECTRONS:
                        values = line.split()
                        id = int(values[0])
                        type_coord_etc = [int(values[1])]
                        for v in map(float, values[2:]):
                            type_coord_etc.append(v)
                        self._handler.process_atoms(id, type_coord_etc)
                    elif state is _ReadState.IONS:
                        values = line.split()
                        self._handler.process_velocities(
                            int(values[0]),
                            map(float, values[1:]))
                    elif state is _ReadState.CELL:
                        values = line.split()
                    elif state is _ReadState.ATOMIC_SPECIES:
                        values = line.split()
                        self._handler.process_velocities(
                            int(values[0]),
                            map(float, values[1:]))
                    elif state is _ReadState.K_POINTS:
                        values = line.split()
                        self._handler.process_velocities(
                            int(values[0]),
                            map(float, values[1:]))
                    elif state is _ReadState.ATOMIC_POSITIONS:
                        values = line.split()

                        self._handler.process_number_atom_types(
                            number_types)
                        self._handler.process_velocities(
                            int(values[0]),
                            map(float, values[1:]))
                        self._handler.process_masses(
                            int(values[0]),
                            int(values[1]))

                    else:
                        continue
            except Exception:
                print("problem with line number=", line_number, line)
                raise
        self._handler.end()


class _ReadState(Enum):
    UNKNOWN, UNSUPPORTED, \
	CONTROL_BEGIN, CONTROL, \
	SYSTEM_BEGIN, SYSTEM, \
	ELECTRONS_BEGIN, ELECTRONS, \
	IONS_BEGIN, IONS,	\
	CELL_BEGIN, CELL, \	
        ATOM_TYPES_BEGIN, ATOM_TYPES, \  #called 'ATOMIC_SPECIES' in qe input files
        ATOMS, \
        ATOMS_BEGIN = range(9),\  #this is for atomic positions maybe? check in cuds
	K_POINTS_BEGIN, K_POINTS, \
	ATOMIC_POSITIONS_BEGIN, ATOMIC_POSITIONS

    @staticmethod
    def get_state(current_state, line):
        """ Reads line and returns state
        """

        new_state = current_state

 	#if we are at state_BEGIN, new state is state not state_BEGIN
 	#if we are at state, new state is state
        if current_state is _ReadState.CONTROL_BEGIN:
            new_state = _ReadState.CONTROL
        elif current_state is _ReadState.SYSTEM_BEGIN:
            new_state = _ReadState.SYSTEM
        elif current_state is _ReadState.ELECTRONS_BEGIN: 
            new_state = _ReadState.ELECTRONS
        elif current_state is _ReadState.IONS_BEGIN:
            new_state = _ReadState.IONS
        elif current_state is _ReadState.CELL_BEGIN: 
            new_state = _ReadState.CELL
        elif current_state is _ReadState.ATOM_TYPES_BEGIN:
            new_state = _ReadState.ATOM_TYPES
        elif current_state is _ReadState.K_POINTS_BEGIN:
            new_state = _ReadState.K_POINTS
        elif current_state is _ReadState.ATOMIC_POSITIONS_BEGIN:
            new_state = _ReadState.ATOMIC_POSITIONS

        else:
            new_state = _ReadState.UNKNOWN

        if "&control" in line:
            new_state = _ReadState.CONTROL_BEGIN
        elif "&system" in line:
            new_state = _ReadState.SYSTEM_BEGIN
        elif "&electrons" in line:
            new_state = _ReadState.ELECTRONS_BEGIN
        elif "&ions" in line:
            new_state = _ReadState.IONS_BEGIN
        elif "&cell" in line:
            new_state = _ReadState.CELL_BEGIN
        elif "ATOMIC_SPECIES" in line:
            new_state = _ReadState.ATOM_TYPES_BEGIN
        elif "K_POINTS" in line:
            new_state = _ReadState.K_POINTS_BEGIN
        elif "ATOMIC_POSITIONS" in line:
            new_state = _ReadState.ATOMIC_POSITIONS_BEGIN

        elif new_state is _ReadState.UNKNOWN:
            new_state = _ReadState.UNSUPPORTED

        return new_state
