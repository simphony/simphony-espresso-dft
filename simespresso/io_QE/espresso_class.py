from enum import Enum
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.particles import Particle, Particles
from simphony.cuds.lattice import Lattice
import numpy as np
import subprocess
#from simphony.cuds.abstractparticles import ABCParticles

import uuid
import logging

import os.path

logging.basicConfig(level=logging.DEBUG)



class qe_functions(object):
    '''
    functions for reading and writing quantum espresso input and output files
    '''

    def __init__(self):
        self.SP = DataContainer()   #System Parameters and Conditions
#        SP = self.data_container
        self.pc =  Particles('quantum_espresso_particles')



    def ReadEspressoOutputFile(self,file_name):
        '''
        This function parses  Espresso output files which usually will have a name like 'name.charge'.
        This file has the structure shown at http://phycomp.technion.ac.il/~sbgrosso/QE_charge_density/node18.html
        The first line has 8 numbers, the three first numbers are the size of the grid in the three spatial directions,
        the three following are the same repeated.
        The seventh number corresponds to the number of atoms and the last one to the number of different type of atoms.
        Second line - the first number is the type of Bravais lattice,
        and the three following numbers are the celldimensions defined in the input file of pw.x.
        The charge density is defined for each point of the grid starting after the atom definitions.
        :param file_name: name of the espresso output file
        :return:
        '''

        if not(os.path.exists(file_name)):
            logging.debug("file "+str(file_name)+" not found")
            return(1)

        with open(file_name, 'r') as f:
            file_iter = iter(f)
            #first line blank
            line = file_iter.next()

            #read first four lines of the header
            try:
                    #2nd line : gridsize x,y,z twice , natoms natomtypes
                    line = file_iter.next()
                    logging.debug('read line2:'+str(line))
                    values = line.split()
                    ints = [int(val) for val in values]
                    n_lattice_points = ints[0:3]
                    n_atoms = ints[6]
                    print('nlattice points:'+str(n_lattice_points)+' n_atoms:'+str(n_atoms))

                    #3rd line : bravais lattice, celldm[0],[1],[2]
                    line = file_iter.next()
                    logging.debug('read line3:'+str(line))
                    values = line.split()
                    floats = [float(val) for val in values]
                    bravais = int(floats[0])
                    celldm = floats[1:4]
    #                SP[CUBA.ORIGINAL_POSITION] = [celldm[0],celldm[1],celldm[2]]
                    print('bravais:'+str(bravais)+' celldm:'+str(celldm))
                    # see http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#idp82064
                    if bravais == 0:
                        Ltype = 'free'
                    if bravais == 1:
                        Ltype = 'cubic P(sc)'
                    if bravais == 2:
                        Ltype = 'cubic F(fcc)'
                    if bravais == 3:
                        Ltype = 'cubic I(bcc)'
                    if bravais == 4:
                        Ltype = 'orthorhombic'
                    if bravais == 5:
                        Ltype = 'trigonal'
                    if bravais == 6:
                        Ltype = 'tetragonal P(st)'
                    if bravais == 7:
                        Ltype = 'tetragonal I(bct)'
                    if bravais == 8:
                        Ltype = 'orthorhombic P'
                    if bravais == 9:
                        Ltype = 'orthorhombic bco'
                    if bravais == 10:
                        Ltype = 'orthorhombic fc'
                    if bravais == 11:
                        Ltype = 'orthorhombic bc'
                    if bravais == 12:
                        Ltype = 'monoclinic P unique axis c'
                    if bravais == 13:
                        Ltype = 'monoclinic base centered'
                    if bravais == 14:
                        Ltype = 'triclinic'

                    self.L = Lattice('quantum espresso lattice',Ltype,celldm,n_lattice_points,[0,0,0])

                    #4th line - don't care
                    line = file_iter.next()
                    logging.debug('read line4:'+str(line))
                    values = line.split()
                    floats = [float(val) for val in values]

            except StopIteration:
                print('eof reached')
                logging.warning('eof reached')
            except :
                #print("problem with line", line)
                logging.warning("problem with line", line)
                return
            #TODO skip all the atom definition lines, eg look for alphabetic characters at beginning of line
            logging.debug('skipping '+str(n_atoms)+' lines')
            for i in range(0,n_atoms+1):
                    line = file_iter.next()

            logging.debug('read '+str(n_lattice_points)+' lattice point lines')
            self.read_densities(n_lattice_points,file_iter,aviz_filename='avizout.xyz')
            #put pc into dc!



    def running_index_to_node_index(self,index,n_latticepoints):
        node_z = index/(n_latticepoints[0]*n_latticepoints[1])
        node_y = (index/n_latticepoints[0])%n_latticepoints[1]
        node_x = index%n_latticepoints[0]
        return [node_x,node_y,node_z]

    def read_densities(self,n_latticepoints,file_iter,aviz_filename=False):
        charge_density = self.read_xyz(n_latticepoints,file_iter)
        charge_as_list = charge_density.flatten()
        nodelist = []
        for charge_index, charge in enumerate(charge_as_list):
            node_index = self.running_index_to_node_index(charge_index,n_latticepoints)# TODO convert charge_index to node index
            node = self.L.get_node(node_index)
            node.data[CUBA.MASS] = charge
            nodelist.append(node)
        self.L.update_nodes(nodelist)  #we have to ask the lattice to update the changed node


        if aviz_filename:
            self.write_aviz_output(charge_density,aviz_filename)
#            self.write_aviz_output(charge_density,aviz_filename,L.base_vect)
        #return L

    #    iterator = L.iter_nodes(n_latticepoints)
    #    for i in range(0,n_latticepoints[0]):
    #        for j in range(0,n_latticepoints[1]):
    #            for k in range(0,n_latticepoints[2]):
    #                L.data = charge_density[i,j,k]
     #               iterator.next()
    def write_aviz_output(self,xyz_array,aviz_xyzfile):
        base_vector = self.L.base_vect
        n_elements = np.shape(xyz_array)
        n_size = np.size(xyz_array)

        with open(aviz_xyzfile, 'w') as f:
            try:
                line = str(n_size)+'\n'  #calculation
                f.write(line)
                line = 'simphony to aviz xyz file\n'  #calculation
                f.write(line)
                for i in range(0,n_elements[0]):
                    x=i*base_vector[0]
                    for j in range(0,n_elements[1]):
                        y=j*base_vector[0]
                        for k in range(0,n_elements[2]):
                            z=k*base_vector[0]
                            density = xyz_array[i,j,k]
                            line = 'C '+str(x)+' ' +str(y)+' '+str(z)+' '+str(density)+'\n'
                            f.write(line)
            except:
                print('error writing aviz file')
        return

    def read_xyz(self,n_latticepoints,file_iter):
        line_numer = 0
        line = file_iter.next()
        x_points = n_latticepoints[0]
        y_points = n_latticepoints[1]
        z_points = n_latticepoints[2]
        x_count = 0
        y_count = 0
        z_count = 0
        charge_density = np.zeros([x_points,y_points,z_points])
        line_number = 0
        try:
            while line is not None:
                line_number = line_number + 1
                print('x{0} y{1} z{2} line:{3}'.format(x_count,y_count,z_count,str(line)))
                values = line.split()
                charges = [float(val) for val in values]
    #            print('charges:'+str(charges))

                for charge in charges:
                    charge_density[x_count,y_count,z_count] = charge
                    x_count+=1
                    if x_count== x_points:
                        x_count = 0
                        y_count+=1
                        if y_count == y_points:
                            y_count = 0
                            z_count +=1
                            if z_count == z_points:
                                break

                    line = file_iter.next()

        except StopIteration:
            print('EOF')
        if x_count<x_points or y_count<y_points or z_count<z_points:
            logging.debug('Got fewer points than expected:read {0} of {1} x, {2} of {3} y, {4} of {5} z'.format(x_count,x_points,y_count,y_points,z_count,z_points))
        return charge_density


    def WriteEspressoInputFile(self,file_name):
        """
        :param file_name: name of the input file to write
        :return:
        """
        SP = self.SP
        pc = self.pc
        #write parameters for a particular working input file

        print('attempting to write '+file_name)
        with open(file_name, 'w') as f:
            try:
                #CONTROL section
                #apparently a comma is not required at the end of every line
                line = '&CONTROL\n'  #calculation
                f.write(line)
                if CUBA.TORQUE in SP:
                    line = '\t calculation=\''+str(SP[CUBA.TORQUE])+'\'\n'   #calculation
                    f.write(line)
                if CUBA.ZETA_POTENTIAL in SP:
                    line = '\t restart_mode=\''+str(SP[CUBA.ZETA_POTENTIAL])+'\'\n' #restart_mode
                    f.write(line)
                if CUBA.YOUNG_MODULUS in SP:
                    line = '\t pseudo_dir=\''+str(SP[CUBA.YOUNG_MODULUS])+'\'\n'  #pseudo dir
                    f.write(line)
                if CUBA.VOLUME_FRACTION in SP:
                    line = '\t prefix=\''+str(SP[CUBA.VOLUME_FRACTION])+'\'\n'  #prefix
                    f.write(line)
                if CUBA.AMPHIPHILICITY in SP:
                    line = '\t tprnfor='+str(SP[CUBA.AMPHIPHILICITY])+'\n'  #tprnfor
                    f.write(line)
                if CUBA.NUMBER_OF_TIME_STEPS in SP:
                    line = '\t max_seconds='+str(int(SP[CUBA.NUMBER_OF_TIME_STEPS]))+'\n'
                    f.write(line)
                if CUBA.DIRECTION in SP:
                    line = '\t outdir=\''+str(SP[CUBA.DIRECTION])+'\'\n'  #outdir
                    f.write(line)
                line = '/\n'
                f.write(line)

                #SYSTEM section
                line = '&SYSTEM\n'  #calculation
                f.write(line)
                if CUBA.ROLLING_FRICTION in SP:
                    line = '\t ibrav='+str(SP[CUBA.ROLLING_FRICTION])+'\n'  #outdir
                    f.write(line)
                if CUBA.ORIGINAL_POSITION in SP:
                    line = '\t celldm(1)='+str(SP[CUBA.ORIGINAL_POSITION][0])+'\n'
                    f.write(line)
                    line = '\t celldm(2)='+str(SP[CUBA.ORIGINAL_POSITION][1])+'\n'
                    f.write(line)
                    line = '\t celldm(3)='+str(SP[CUBA.ORIGINAL_POSITION][2])+'\n'
                    f.write(line)

                 # here goes nat and ntype
                n_atoms = self.count_particles()
                line = '\t nat='+str(n_atoms)+'\n'  #outdir
                f.write(line)
                if CUBA.SCALING_COEFFICIENT in SP:
                    line = '\t ntyp='+str(SP[CUBA.SCALING_COEFFICIENT])+'\n'  #outdir
                    f.write(line)



                if CUBA.LN_OF_RESTITUTION_COEFFICIENT in SP:
                    line = '\t ecutwfc='+str(SP[CUBA.LN_OF_RESTITUTION_COEFFICIENT])+'\n'
                    f.write(line)
                if CUBA.POISSON_RATIO in SP:
                    line = '\t ecutrho='+str(SP[CUBA.POISSON_RATIO])+'\n'
                    f.write(line)
                if CUBA.LATTICE_SPACING in SP:
                    line = '\t input_dft=\''+str(SP[CUBA.LATTICE_SPACING])+'\'\n'  #outdir
                    f.write(line)
                line = '/\n'
                f.write(line)


                #ELECTRONS section
                line = '&ELECTRONS\n'  #calculation
                f.write(line)
                if CUBA.SMOOTHING_LENGTH in SP:
                    line = '\t mixing_mode=\''+str(SP[CUBA.SMOOTHING_LENGTH])+'\'\n'
                    f.write(line)
                if CUBA.PHASE_INTERACTION_STRENGTH in SP:
                    line = '\t mixing_beta='+str(SP[CUBA.PHASE_INTERACTION_STRENGTH])+'\n'
                    f.write(line)
                if CUBA.DEBYE_LENGTH in SP:
                    line = '\t conv_thr='+str(SP[CUBA.DEBYE_LENGTH])+'\n'
                    f.write(line)
                line = '/\n'
                f.write(line)

                #IONS section
                line = '&IONS\n'  #calculation
                f.write(line)
                line = '/\n'
                f.write(line)

                #CELL section
                line = '&CELL\n'  #calculation
                f.write(line)
                line = '/\n'
                f.write(line)

                #CELL section
                line = 'ATOMIC_SPECIES\n'  #calculation
                f.write(line)

                #ATOMIC SPECIES
                if CUBA.CHEMICAL_SPECIE in SP:
                    for i in range(0,len(SP[CUBA.CHEMICAL_SPECIE])):
                        line = '\t'+str(SP[CUBA.CHEMICAL_SPECIE][i])+' '+\
                               str(SP[CUBA.MASS][i])+' '+str(SP[CUBA.FRICTION_COEFFICIENT][i])+'\n'
                        f.write(line)
                line = '\n'
                f.write(line)

                #K POINTS
                if CUBA.PROBABILITY_COEFFICIENT in SP:
                    line ='K_POINTS '+str(SP[CUBA.PROBABILITY_COEFFICIENT])+'\n'
                    f.write(line)
                if CUBA.EQUATION_OF_STATE_COEFFICIENT in SP:
                    line = ''
                    for k_point in SP[CUBA.EQUATION_OF_STATE_COEFFICIENT]:
                        line = line+str(k_point)+' '
                    f.write(line)
                f.write('\n\n')

                if CUBA.KINEMATIC_VISCOSITY in SP:
                    line ='ATOMIC_POSITIONS '+str(SP[CUBA.KINEMATIC_VISCOSITY])+'\n'
                    f.write(line)


                multiplier = 10**10  #QE wants coords. in Angstroms
                for particle in pc.iter_particles():
                    atom_type = particle.data
                 #   specie = atom_type[CUBA.CHEMICAL_SPECIE]
                    atom = atom_type[CUBA.CHEMICAL_SPECIE]
                    print('atom:'+str(atom)+' data:'+str(atom_type)+' ')
                    line =str(atom)+' '+str(multiplier*particle.coordinates[0])+' '+\
                          str(multiplier*particle.coordinates[1])+' '+str(multiplier*particle.coordinates[2]) + '\n'
                    f.write(line)

            except:
                    ('error in write block of WriteEspressoInputFile')
                    raise
        print('finished writing file')
        f.closed

    def kluge(self,f):
        line =' C 1.0 2.0 3.0 \n'
        f.write(line)
        line =' C 3.0 3.0 4.0 \n'
        f.write(line)
        line =' C 3.0 4.0 5.0 \n'
        f.write(line)
        line =' C 4.0 5.0 6.0 \n'
        f.write(line)
        return

    def WriteEspressoPPFile(self,ppfilename="testpp.in"):
        '''
        this writes an auxiliary required file determined the plot parameters
        :return:
        '''
        outdir = './'
        plotfile = 'output.charge'
        outfile = 'density.dat'
        lines=['&inputpp',
               'prefix =\'qe_output\'',
               'filplot=\''+plotfile+'\'',
               'plot_num=0',
               'outdir = '+outdir,
               '/',
                '&plot',
               'nfile = 1',
                'filepp(1) = \''+plotfile+'\'',
               'weight(1) = 1.0',
                'iflag = 3',
                'output_format = 6',
                'fileout = \''+outfile+'\''  ]
        with open(ppfilename,'w') as pp:
            for line in lines:
                print'line:'+str(line)
                pp.write(str(line)+'\n')
        print('finished writing file')

    def ReadEspressoInputFile(self,file_name):
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
        self.celldm=[None,None,None]
        state = _ReadState.UNKNOWN
    #    BC = DataContainer()   #boundary conditions
  #      CM = DataContainer()   #computational method
        SP = self.SP   #System Parameters and Conditions
  #      SD = DataContainer()  #state data
        pc = self.pc
 #       dc = DataContainer()
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
                    logging.debug('read line:'+str(line))
                    # skip blank lines
     #               if not line.strip():
     #                   continue

                    state = _ReadState.get_state(state,line)
                    if state is _ReadState.CONTROL:
                        print('reading control section')
                        line = self.process_control(file_iter,SP)
                        continue
                    elif state is _ReadState.SYSTEM:
                        print('reading system')
                        line = self.process_system(file_iter,SP)
                        continue
                    elif state is _ReadState.ELECTRONS:
                        print('reading electrons')
                        line = self.process_electrons(file_iter,SP)
                        continue
                    elif state is _ReadState.IONS:
                        print('reading ions')
                        line = self.process_ions(file_iter,SP)
                        continue
                    elif state is _ReadState.CELL:
                        print('reading ions')
                        line = self.process_ions(file_iter,SP)
                        continue
                    elif state is _ReadState.ATOMIC_SPECIES:
                        print('reading atomic species')
                        line = self.process_atomic_species(file_iter,SP)
                        continue
                    elif state is _ReadState.K_POINTS:
                        print('reading k points')
                        values = line.split()
                        line = self.process_k_points(file_iter,SP,mode=values[1])
                        continue
                    elif state is _ReadState.ATOMIC_POSITIONS:
                        print('reading atomic positions')
                        values = line.split()
                        pc = self.process_atomic_positions(file_iter,pc,SP,units=values[1])

                    break

                    line = file_iter.next()
            except StopIteration:
                print('eof reached')
            except Exception:
                print("problem with line number=", line_number, line)
                return
            #put pc into dc!

        return pc


    def process_control(self,f,SP):
        print('processing control section')
        line = f.next()
        while _ReadState.get_state(_ReadState.CONTROL,line) == _ReadState.CONTROL:
            values = [x.strip() for x in line.split('=')]
            logging.debug('line in control section:'+str(line))
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

    def process_system(self,f,SP):
        self.celldm=[None,None,None]
        print('processing system section')
        line = f.next()
        while _ReadState.get_state(_ReadState.SYSTEM,line) == _ReadState.SYSTEM:
            values = [x.strip() for x in line.split('=')]
            logging.debug('line in control section:'+str(line))
            if "ibrav" in line:  #bravais lattice index
                ibrav = int(values[1])
                SP[CUBA.ROLLING_FRICTION] = ibrav
            elif "celldm(1)" in line:
                self.celldm[0] = float(values[1])
    #            SP[CUBA.ORIGINAL_POSITION] = celldm[0]
            elif "celldm(2)" in line:
                self.celldm[1] = float(values[1])
    #            SP[CUBA.ORIGINAL_POSITION][1] = celldm[1]
            elif "celldm(3)" in line:
                self.celldm[2] = float(values[1])
                SP[CUBA.ORIGINAL_POSITION] = [self.celldm[0],self.celldm[1],self.celldm[2]]
            #    SP[CUBA.LATTICE_VECTORS] = celldm

            elif "nat" in line:
                n_atoms = int(values[1])
            elif "ntyp" in line:
                n_atom_types = int(values[1])
                SP[CUBA.SCALING_COEFFICIENT] = n_atom_types
            elif "ecutwfc" in line:
                ecutwfc = float(values[1])
                SP[CUBA.LN_OF_RESTITUTION_COEFFICIENT] = ecutwfc
            elif "ecutrho" in line:
                ecutrho = float(values[1]) #maybe int
                SP[CUBA.POISSON_RATIO] = ecutrho
            elif "input_dft" in line:
                input_dft = values[1]
                SP[CUBA.LATTICE_SPACING] = input_dft
            line = f.next()
        return line


    def process_electrons(self,f,SP):
        print('processing eletrons section')
        line = f.next()
        while _ReadState.get_state(_ReadState.ELECTRONS,line) == _ReadState.ELECTRONS:
            values = [x.strip() for x in line.split('=')]
            logging.debug('line in electrons section:'+str(line))
            if "mixing_mode" in line:
                mixing_mode = values[1]
                SP[CUBA.SMOOTHING_LENGTH] = mixing_mode
            elif "mixing_beta" in line:
                mixing_beta = float(values[1])
                SP[CUBA.PHASE_INTERACTION_STRENGTH] = mixing_beta
            elif "conv_thr" in line:
                conv_thr = values[1] #numbers like 1.0d-7 might have to be converted to float
                SP[CUBA.DEBYE_LENGTH] = conv_thr
            line = f.next()
        return line

    def process_ions(self,f,SP):
        print('processing ions section')
        line = f.next()
        line = f.next()
        return line

    def process_cell(self,f,SP):
        print('processing cell section')
        line = f.next()
        line = f.next()
        return line

    def process_atomic_species(self,f,SP):
        print('processing atomic species section')
        line = f.next()
        SP[CUBA.CHEMICAL_SPECIE]=[]
        SP[CUBA.MASS]=[]
        SP[CUBA.FRICTION_COEFFICIENT]=[]
        while _ReadState.get_state(_ReadState.ATOMIC_SPECIES,line) == _ReadState.ATOMIC_SPECIES:
            values = line.split()
            logging.debug('line in atomic species section:'+str(line))
    #            print('atomtypes:'+str(self.atomtypes))

            if len(values)>0:
                if values[0] in atomtypes:
                    print("atom type:"+values[0])
                    #self.dc(CHEMICAL_SPECIE = values[0])
                    SP[CUBA.CHEMICAL_SPECIE].append(values[0])
                    mass = float(values[1])
                    SP[CUBA.MASS].append(mass)
                    potential_file = values[2]
                    SP[CUBA.FRICTION_COEFFICIENT].append(potential_file)
            line = f.next()
        return line

    def process_k_points(self,f,SP,mode='automatic'):
        #skip line
        print('processing k_points section')
        line = f.next()
    #    SP[CUBAExtension.K_POINTS_MODE] = mode
        i = 0

        SP[CUBA.PROBABILITY_COEFFICIENT] = mode

        while _ReadState.get_state(_ReadState.K_POINTS,line) == _ReadState.K_POINTS:
    #        print('line:'+str(line))
            values = line.split()
            if len(values):
                K_points = (values)
                print('k points:'+str(K_points))
                SP[CUBA.EQUATION_OF_STATE_COEFFICIENT] = K_points

            line = f.next()
        return line


    #I am not sure whether to use datacontainer or particlecontainer - maybe PC goes in DC?
    def process_atomic_positions(self,f,pc,SP,units='(angstrom)'):
        print('processing atomic_positions section')
        try:
            line = f.next()
            SP[CUBA.KINEMATIC_VISCOSITY] = units
            particle_list = []
            while _ReadState.get_state(_ReadState.ATOMIC_POSITIONS,line) == _ReadState.ATOMIC_POSITIONS:
                logging.debug('line in atomic positions section:'+str(line))
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
                    print('uid:'+str(p.uid))
                    p.data[CUBA.CHEMICAL_SPECIE] = atomtype
                    particle_list.append(p)
                    try:
                        line = f.next()
                    except StopIteration:
                        print('EOF')
                        break
            if len(particle_list):
                pc.add_particles(particle_list)
                  #  part_container.add_particle(p)
                n=self.count_particles()
                print('n_particles='+str(n))
        except StopIteration:
            return('EOF')

        return pc

    def count_particles(self):
        n=0
        for particle in self.pc.iter_particles():
            n+=1
        return n

    def start_qe(self,name_in,name_out,path_to_espresso='/usr/bin/pw.x',mpi=False,mpi_Nprocessors=2):
#        name_in = './test_pw.in'
#        name_out = './test_pw.out'
        if path_to_espresso is None:
            path_to_espresso = '/usr/bin/pw.x'   #this appears to be the default install location for espresso
        if mpi:
            command = 'mpirun -np '+str(mpi_Nprocessors)+' '+path_to_espresso+' < '+name_in +' > '+name_out
        else:
            command = path_to_espresso+' < '+name_in +' > '+name_out

       # command = '/usr/bin/pw.x < '+name_in+' > '+name_out
        print('qe wrapper attempting to run: '+command)
        #alternative would be to use subprocess.check_call()  -  however this would give the same info as the
        # try/except, while taking twice as long in the case of success, iiuc
        try:
            subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()
        except:
            e = sys.exc_info()[0]
            print( "<p>Error: %s</p>" % e)

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
            new_state = _ReadState.IONS
        elif "&CELL" in line:
            new_state = _ReadState.CELL
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
    wrapper = qe_wrapper()
#    filename = 'xyzoutput.txt.bak'
    filename = '../../examples/input_pw.in'
    print('started parsing qe input file '+str(filename))
    wrapper.ReadEspressoInputFile(filename)

    new_inputfilename =  '../../examples/new_input_pw.in'
    wrapper.WriteEspressoInputFile(new_inputfilename)

    ppfilename = 'testpp.in'
    print('started writing qe input file '+str(ppfilename))
    wrapper.WriteEspressoPPFile()

    filename = 'xyzoutput.txt'
    print('started parsing file '+str(filename))
    wrapper.ReadEspressoOutputFile(filename)


