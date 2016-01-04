import logging
import os
import os.path
import subprocess
import sys
import numpy as np
from enum import Enum
import math

from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_particles import ABCParticles
from simphony.cuds.lattice import Lattice
from simphony.cuds.particles import Particle, Particles

from qeCubaExtensions import qeCUBAExtension
import data_handler

logging.basicConfig(level=logging.DEBUG)


class QeWrapper(ABCModelingEngine):
    '''
    functions for reading and writing quantum espresso input and output files
    wrapper must define:
        add_dataset(self, container):
        remove_dataset(self, name):
        get_dataset(self, name):
        get_dataset_names(self):
        iter_datasets(self, names=None):
        run(self)
        _combine(data_container, data_container_extension)
           ?? not sure about last one - it doe not appear in
           common/abc_modeling_engine
    '''

    def __init__(self):
        self.datahandler=data_handler.qe_data_handler()


    def run(self):
        print('starting qe engine')
        self.write_espresso_input_file(self.input_pwname)
        print('path to espresso:'+self.path_to_espresso)
#        print(which(self.path_to_espresso))
        if not which(self.path_to_espresso):
            logging.debug('no path to espresso')
            raise ValueError(
                'espresso command not found (looking for '
                            + self.path_to_espresso+')')
        if self.mpi:
            command = 'mpirun -np ' + str(self.mpi_Nprocessors) + ' ' + \
                      self.path_to_espresso + ' < ' + self.pwname_in + ' > ' \
                      + self.output_filename
        else:
            command = self.path_to_espresso + ' < ' + self.input_pwname + ' > ' \
                      + self.output_filename
        logging.debug('attempting to run command: ' + command)
        try:
            subprocess.check_call(command, shell=True,
                stdout=subprocess.PIPE).stdout.read()
#            subprocess.Popen(command, shell=True,
#                             stdout=subprocess.PIPE).stdout.read()
        except:
            e = sys.exc_info()[0]
            logging.debug('espresso command gave error %s' % e)
            raise ValueError('espresso command gave error %s' % e)

    def add_dataset(self,container):
        """Add a CUDS container
        Parameters
        ----------
        container : {ABCParticles}
            The CUDS container to add to the engine.
        Raises
        ------
        TypeError:
            If the container type is not supported (i.e. ABCLattice, ABCMesh).
        ValueError:
            If there is already a dataset with the given name.
        """
        if not isinstance(container, ABCParticles):
            #if we allow lattice or mesh this allows for charge density
            #also would need several, one for each energy state
            # currently lattice lacks  enough metadata for this
            raise TypeError(
                "This type of dataset container is not supported")

        if container.name in self.get_dataset_names():
            raise ValueError(
                'Particle container \'{}\' already exists'.format(
                    container.name))
        else:
            self.datasets[container.name] = container

    def remove_dataset(self,dataset_name):
        """remove dataset
        Parameters
        ----------
        dataset_name : name container to remove
        Raises
        ------
        ValueError:
            If there is no dataset with the given name.
        """
        if not dataset.name in self.dataset_names:
            raise ValueError(
                'Particle container \'{}\' does not exist'.format(
                    dataset_name))
        else:
            del self.datasets[dataset.name]

    def get_dataset(self,dataset_name):
        """Get a CUDS container
        Parameters
        ----------
        container_name : string
            Name of the CUDS container to get.
        Raises
        ------
        ValueError:
            If there is no dataset with the given name.
        """
        if not dataset_name in self.get_dataset_names():
            logging.debug('Container name '+str(dataset_name)+' not found.')
            return None
        else:
            dataset = self.datasets[dataset_name]
            return dataset

    def get_dataset_names(self):
        """Get names of datasets in container
        """
        dataset_names = [name for name in self.datasets]
        return dataset_names

    def iter_datasets(self, names=None):
        """ Returns an iterator over a subset or all of the containers.
        Parameters
        ----------
        names : sequence of str, optional
            names of specific containers to be iterated over. If names is not
            given, then all containers will be iterated over.
        """
        if names is None:
            for name in self.get_dataset_names:
                yield self.get_dataset(name)
        else:
            for name in names:
                if name in self.get_dataset_names:
                    yield self.get_dataset(name)
                else:
                    raise ValueError(
                        'Particle container \'{}\` does not exist'.format(
                            name))


    def read_espresso_output_file(self, file_name):
        '''
        This function parses  Espresso output files which usually will have
        a name like 'name.charge'.
        This file has the structure shown at
        http://phycomp.technion.ac.il/~sbgrosso/QE_charge_density/node18.html
        The first line has 8 numbers, the three first numbers are the size of
        the grid in the three spatial directions, the three following are the
        same repeated.
        The seventh number corresponds to the number of atoms and the last
        one to the number of different type of atoms.
        Second line - the first number is the type of Bravais lattice,
        and the three following numbers are the celldimensions defined in the
        input file of pw.x.
        The charge density is defined for each point of the grid starting
        after the atom definitions.
        :param file_name: name of the espresso output file
        :return:
        '''
        if not (os.path.exists(file_name)):
            logging.debug("file " + str(file_name) + " not found")
            return (1)

        with open(file_name, 'r') as f:
            file_iter = iter(f)

            try:
                # read first four lines of the header - first line blank
                line = file_iter.next()

                # 2nd line : gridsize x,y,z twice , natoms natomtypes
                line = file_iter.next()
                logging.debug('read line2:' + str(line))
                values = line.split()
                ints = [int(val) for val in values]
                for val in ints:
                    if not isinstance(val, int):
                        logging.debug('got noninteger value in output, line 2')
                        return None
                n_lattice_points = ints[0:3]
                n_atoms = ints[6]
                logging.debug(
                    'n_points:' + str(n_lattice_points) +
                    ' n_atoms:' + str(n_atoms))

                # 3rd line : bravais lattice, celldm[0],[1],[2]
                line = file_iter.next()
                logging.debug('read line3:' + str(line))
                values = line.split()
                floats = [float(val) for val in values]
                for val in floats:
                    if not isinstance(val, float) and not \
                            isinstance(val, int):
                        logging.debug('got non-float/int in input line 3')
                        return None
                bravais = int(floats[0])
                celldm = floats[1:4]
                logging.debug(
                    'bravais:' + str(bravais) +
                    ' celldm:' + str(celldm))
                # see http://www.quantum-espresso.org/
                # wp-content/uploads/Doc/INPUT_PW.html#idp82064
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

                self.L = Lattice('quantum espresso lattice',
                                 Ltype, celldm, n_lattice_points, [0, 0, 0])
                self.BC.lattice = self.L
                # 4th line - don't care
                line = file_iter.next()
                logging.debug('read line4:' + str(line))
            except StopIteration:
                print('eof reached')
                logging.warning('eof reached')
            except:
                # print("problem with line", line)
                if (line):
                    print(str(line))
                    logging.warning("problem with line", str(line))
                else:
                    logging.warning("no line obtained")
                return
            # TODO skip all the atom definition lines,
            # eg look for alphabetic characters at beginning of line
            logging.debug('skipping ' + str(n_atoms) + ' lines')
            for i in range(0, n_atoms + 1):
                line = file_iter.next()

            logging.debug('read ' + str(n_lattice_points) +
                          ' lattice point lines')
            self.read_densities(n_lattice_points, file_iter,
                                aviz_filename='avizout.xyz')

    def running_index_to_node_index(self, index, n_latticepoints):
        node_z = index / (n_latticepoints[0] * n_latticepoints[1])
        node_y = (index / n_latticepoints[0]) % n_latticepoints[1]
        node_x = index % n_latticepoints[0]
        return [node_x, node_y, node_z]

    def read_densities(self, n_latticepoints, file_iter, aviz_filename=False):
        charge_density = self.read_xyz(n_latticepoints, file_iter)
        charge_as_list = charge_density.flatten()
        nodelist = []
        for charge_index, charge in enumerate(charge_as_list):
            node_index = self.running_index_to_node_index(charge_index,
                                                          n_latticepoints)
            # TODO convert charge_index to node index
            node = self.L.get_node(node_index)
            # TODO check if this is ok
            node.data[CUBA.MASS] = charge
            nodelist.append(node)
        self.L.update_nodes(nodelist)
        # we have to ask the lattice to update the changed node

        if aviz_filename:
            self.write_aviz_output(charge_density, aviz_filename)

    def write_aviz_output(self, xyz_array, aviz_xyzfile):
        base_vector = self.L.base_vect
        n_elements = np.shape(xyz_array)
        n_size = np.size(xyz_array)

        with open(aviz_xyzfile, 'w') as f:
            try:
                line = str(n_size) + '\n'
                f.write(line)
                line = 'simphony to aviz xyz file\n'
                f.write(line)
                for i in range(0, n_elements[0]):
                    x = i * base_vector[0]
                    for j in range(0, n_elements[1]):
                        y = j * base_vector[0]
                        for k in range(0, n_elements[2]):
                            z = k * base_vector[0]
                            density = xyz_array[i, j, k]
                            line = 'C ' + str(x) + ' ' + str(y) + ' ' + str(z) + \
                                   ' ' + str(density) + '\n'
                            f.write(line)
            except:
                print('error writing aviz file')
        return

    def read_xyz(self, n_latticepoints, file_iter):
        line = file_iter.next()
        x_points = n_latticepoints[0]
        y_points = n_latticepoints[1]
        z_points = n_latticepoints[2]
        x_count = 0
        y_count = 0
        z_count = 0
        charge_density = np.zeros([x_points, y_points, z_points])
        line_number = 0
        try:
            while line is not None:
                line_number = line_number + 1
                print('x{0} y{1} z{2} line:{3}'
                      .format(x_count, y_count, z_count, str(line)))
                values = line.split()
                charges = [float(val) for val in values]
                #            print('charges:'+str(charges))

                for charge in charges:
                    charge_density[x_count, y_count, z_count] = charge
                    x_count += 1
                    if x_count == x_points:
                        x_count = 0
                        y_count += 1
                        if y_count == y_points:
                            y_count = 0
                            z_count += 1
                            if z_count == z_points:
                                break

                    line = file_iter.next()

        except StopIteration:
            print('EOF')
        if x_count < x_points or y_count < y_points or z_count < z_points:
            logging.debug('Got fewer points than expected:read {0} of {1} x, '
                          '{2} of {3} y, {4} of {5} z'.
                          format(x_count, x_points, y_count,
                                 y_points, z_count, z_points))
        return charge_density

    def write_espresso_input_file(self, file_name,dataset_name = None):
        """
        :param file_name: name of the input file to write
        :return:
        """
        SP = self.SP
        if dataset_name is  None:
            dataset_name = self.get_dataset_names()[0]
        pc = self.get_dataset(dataset_name)

#        SD = self.SD
        # write parameters for a particular working input file

        print('attempting to write ' + file_name)
        try:
            with open(file_name, 'w') as f:
                # CONTROL section
                # apparently a comma is not required at the end of every line
                line = '&CONTROL\n'
                f.write(line)
                if hasattr(self, 'calculation_type'):
                    line = '\t calculation=\'' + str(self.calculation_type) \
                           + '\'\n'
                    f.write(line)
                if hasattr(self, 'restart_mode'):
                    line = '\t restart_mode=\'' + str(self.restart_mode)\
                           + '\'\n'
                    f.write(line)
                if hasattr(self, 'pseudopotential_directory'):
                    line = '\t pseudo_dir=\'' + str(self.pseudopotential_directory) \
                           + '\'\n'
                    f.write(line)
                if hasattr(self, 'pseudopotential_prefix'):
                    line = '\t prefix=\'' + str(self.pseudopotential_prefix) \
                           + '\'\n'
                    f.write(line)
                if hasattr(self, 'tprnfor'):
                    line = '\t tprnfor=' + str(self.tprnfor) + '\n'
                    f.write(line)
                if hasattr(self, 'max_seconds'):
                    line = '\t max_seconds=' + \
                           str(int(self.max_seconds)) + '\n'
                    f.write(line)
                if hasattr(self, 'output_directory'):
                    line = '\t outdir=\'' + str(self.output_directory) + '\'\n'
                    f.write(line)
                line = '/\n'
                f.write(line)

                # SYSTEM section
                line = '&SYSTEM\n'
                f.write(line)
                if hasattr(self, 'ibrav'):
                    line = '\t ibrav=' + str(self.ibrav) + '\n'
                    # outdir
                    f.write(line)
                if hasattr(self, 'celldm'):
                    line = '\t celldm(1)=' + \
                           str(self.celldm[0]) + '\n'
                    f.write(line)
                    line = '\t celldm(2)=' + \
                           str(self.celldm[1]) + '\n'
                    f.write(line)
                    line = '\t celldm(3)=' + \
                           str(self.celldm[2]) + '\n'
                    f.write(line)

                # here goes nat and ntype
                n_atoms = self.count_particles()
                line = '\t nat=' + str(n_atoms) + '\n'
                f.write(line)
                n_atom_types = self.count_atom_types(pc)
                logging.debug('n_atom_types:'+str(n_atom_types))
                line = '\t ntyp=' + str(n_atom_types) + '\n'
                f.write(line)
                if hasattr(self, 'ecutwfc'):
                    line = '\t ecutwfc=' + \
                           str(self.ecutwfc) + '\n'
                    f.write(line)
                if hasattr(self, 'ecutrho'):
                    line = '\t ecutrho=' + str(self.ecutrho) + '\n'
                    f.write(line)
                if hasattr(self, 'input_dft'):
                    line = '\t input_dft=\'' + str(self.input_dft) \
                           + '\'\n'
                    f.write(line)
                line = '/\n'
                f.write(line)

                # ELECTRONS section
                line = '&ELECTRONS\n'
                f.write(line)
                if hasattr(self, 'mixing_mode'):
                    line = '\t mixing_mode=\'' + \
                           str(self.mixing_mode) + '\'\n'
                    f.write(line)
                if hasattr(self, 'mixing_beta'):
                    line = '\t mixing_beta=' + \
                           str(self.mixing_beta) + '\n'
                    f.write(line)
                if hasattr(self, 'convergence_threshold'):
                    logging.debug('conv thresh:'+str(self.convergence_threshold))
                    logdecimal = math.log10(self.convergence_threshold)
                    base = 10 ** (logdecimal - int(logdecimal))
                    logging.debug('exp:'+str(logdecimal))
                    logging.debug('base:'+str(base))
                    line = '\t conv_thr=' + str(base)+'d'+str(int(logdecimal))+ '\n'
                    f.write(line)
                line = '/\n'
                f.write(line)

                # IONS section
                line = '&IONS\n'
                f.write(line)
                line = '/\n'
                f.write(line)

                # CELL section
                line = '&CELL\n'
                f.write(line)
                line = '/\n'
                f.write(line)

                # ATOMIC SPECIES section
                line = 'ATOMIC_SPECIES\n'
                f.write(line)
                # C 12.0107 06-C.GGA.fhi.UPF
                for atomtype in self.what_atom_types():
                    potential_file = self.SP_extension[qeCUBAExtension.PSEUDO_POTENTIAL]
                    #todo - take care of possible isotopes - same atom, different mass
                    particle_mass_list = self.get_particle_masses()
                    if atomtype in particle_mass_list:
                        mass = particle_mass_list[atomtype]
                        logging.debug('mass of {0} is {1}'.format(atomtype,mass))
                    else:
                        logging.warning('could not find mass of particle:'+str(atomtype))
                        continue
                    line =atomtype + ' ' + str(mass)+' ' + potential_file + '\n'
                    logging.debug(line)
                    f.write(line)
                line = '\n'
                f.write(line)

                # K POINTS
                sampling_mode = 'automatic'
#                if hasattr(self,'CM_extension[qeCUBAExtension.K_POINT_SAMPLING_METHOD]'):
                sampling_mode = self.CM_extension[qeCUBAExtension.K_POINT_SAMPLING_METHOD]
                line = 'K_POINTS ' + str(sampling_mode) + '\n'
                f.write(line)

                #figure out how to write this into input.pw
#                line = 'K_POINTS ' + \
 #                      str(self.k_points_mode) + '\n'
  #              f.write(line)

#                if hasattr(self,'CM_extension[qeCUBAExtension.K_POINT_SAMPLING]'):
                kpoints = self.CM_extension[qeCUBAExtension.K_POINT_SAMPLING]
                line = ''
#                if hasattr(self, 'k_point_values'):
                for point in kpoints:
                    line = line + str(point) + ' '
                line = line + '\n'
                logging.debug(line)
                f.write(line)
                line = '\n'
                f.write(line)


                # ATOMIC POSITIONS
                # this will apparently always be angstroms iiuc - jr
                line = 'ATOMIC_POSITIONS ' + '(' + str(self.position_units) + ')'+'\n'
                f.write(line)

                multiplier = 10.0 ** 10
                multiplier = 1.0
                # QE wants coords. in Angstroms
                for particle in pc.iter_particles():
                    atom_type = particle.data
                    #   specie = atom_type[CUBA.CHEMICAL_SPECIE]
                    atom = atom_type[CUBA.CHEMICAL_SPECIE][0]
                    logging.debug('atom:' + str(atom) + ' data:' +
                          str(atom_type) + ' ')
                    line = str(atom) + ' ' + \
                        str(multiplier * particle.coordinates[0]) + ' ' + \
                        str(multiplier * particle.coordinates[1]) + ' ' + \
                        str(multiplier * particle.coordinates[2]) + '\n'
                    f.write(line)
        except:
            ('error in write block of write_espresso_input_file')
            raise
        print('finished writing file')
        f.closed

    def write_espresso_pp_file(self, ppfilename="testpp.in"):
        '''
        this writes an auxiliary required file determined the plot parameters
        :return:
        '''
        outdir = './'
        plotfile = 'output.charge'
        outfile = 'density.dat'
        lines = ['&inputpp',
                 'prefix =\'qe_output\'',
                 'filplot=\'' + plotfile + '\'',
                 'plot_num=0',
                 'outdir = ' + outdir,
                 '/',
                 '&plot',
                 'nfile = 1',
                 'filepp(1) = \'' + plotfile + '\'',
                 'weight(1) = 1.0',
                 'iflag = 3',
                 'output_format = 6',
                 'fileout = \'' + outfile + '\'']
        with open(ppfilename, 'w') as pp:
            for line in lines:
                print'line:' + str(line)
                pp.write(str(line) + '\n')
        print('finished writing file')

    def read_espresso_input_file(self, file_name):
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
        2. convert output (which is a charge density file) into simphony format
        (see charge_density_xyz.cpp)

        Parameters
        ----------
        handler :
           handler will handle the parsed information provided by this class
            calculation_type
            restart_mode
            pseudo_dir
            prefix
            tprnfor
            max_seconds
            outdir
        file_name : name of qe input file
            """
        self.celldm = [None, None, None]
        state = _ReadState.UNKNOWN
        #       BC = DataContainer()   #boundary conditions
        #       CM = DataContainer()   #computational method
        #       SD = DataContainer()  #state data
        pc = self.pc
        #       dc = DataContainer()
        with open(file_name, 'r') as f:
            line_number = 0
            file_iter = iter(f)
            line = file_iter.next()
            try:
                while line is not None:
                    line_number += 1
                    logging.debug('read line:' + str(line))
                    state = _ReadState.get_state(state, line)
                    if state is _ReadState.CONTROL:
                        print('reading control section')
                        line = self.process_control(file_iter)
                        continue
                    elif state is _ReadState.SYSTEM:
                        print('reading system')
                        line = self.process_system(file_iter)
                        continue
                    elif state is _ReadState.ELECTRONS:
                        print('reading electrons')
                        line = self.process_electrons(file_iter)
                        continue
                    elif state is _ReadState.IONS:
                        print('reading ions')
                        line = self.process_ions(file_iter)
                        continue
                    elif state is _ReadState.CELL:
                        print('reading ions')
                        line = self.process_ions(file_iter)
                        continue
                    elif state is _ReadState.ATOMIC_SPECIES:
                        print('reading atomic species')
                        line = self.process_atomic_species(file_iter)
                        continue
                    elif state is _ReadState.K_POINTS:
                        print('reading k points')
                        values = line.split()
                        line = self.process_k_points(file_iter,
                                                     mode=values[1])
                        continue
                    elif state is _ReadState.ATOMIC_POSITIONS:
                        print('reading atomic positions')
                        values = line.split()
                        self.process_atomic_positions(file_iter,
                                                           units=values[1])
                        #  take out pc and use self.PC
                    break

                #                    line = file_iter.next()
            except StopIteration:
                print('eof reached')
            except Exception:
                print("problem with line number=", line_number, line)
                return
        return

    def process_control(self, f):
        print('processing control section')
        line = f.next()
        while _ReadState.get_state(_ReadState.CONTROL, line) == \
                _ReadState.CONTROL:
            values = [x.strip() for x in line.split('=')]
            logging.debug('line in control section:' + str(line))
            if "calculation" in line:
                values = line.split('=')
                calculation_type = values[1]
#               Not sure if this is ok or not...
                self.calculation_type = calculation_type
            elif "restart_mode" in line:
                self.restart_mode = values[1]
# TODO change restart mode to check if qe was interrupted previously
# - should not be a cuba keyword
# set restart mode='restart' if qe was interrupted as described here
# http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#idp27692160
            elif "pseudo_dir" in line:
                pseudo_dir = values[1]
                self.pseudopotential_directory = pseudo_dir
            elif "prefix" in line:
                prefix = values[1]
                self.pseudopotential_prefix = prefix
            elif "tprnfor" in line:
                tprnfor = values[1]
                self.tprnfor = tprnfor
            elif "max_seconds" in line:
                max_seconds = float(values[1])
                self.max_seconds = max_seconds
            elif "outdir" in line:
                outdir = values[1]
                self.output_directory = outdir
            line = f.next()
        return line

    def process_system(self, f):
        self.celldm = [None, None, None]
        print('processing system section')
        line = f.next()
        while _ReadState.get_state(_ReadState.SYSTEM, line) == \
                _ReadState.SYSTEM:
            values = [x.strip() for x in line.split('=')]
            logging.debug('line in control section:' + str(line))
            if "ibrav" in line:
                ibrav = int(values[1])
                self.ibrav = ibrav

            elif "celldm(1)" in line:
                self.celldm = [0, 0, 0]
                self.celldm[0] = float(values[1])
            elif "celldm(2)" in line:
                self.celldm[1] = float(values[1])
            elif "celldm(3)" in line:
                self.celldm[2] = float(values[1])
                self.SP[CUBA.LATTICE_VECTORS] = self.celldm

            elif "nat" in line:
                pass
            elif "ntyp" in line:
                n_atom_types = int(values[1])
                self.n_atom_types = n_atom_types
            elif "ecutwfc" in line:
                ecutwfc = float(values[1])
                self.ecutwfc = ecutwfc
            elif "ecutrho" in line:
                ecutrho = float(values[1])
                # maybe int
                self.ecutrho = ecutrho
            elif "input_dft" in line:
                input_dft = values[1]
                self.input_dft = input_dft
            line = f.next()
        return line

    def process_electrons(self, f):
        print('processing eletrons section')
        line = f.next()
        while _ReadState.get_state(_ReadState.ELECTRONS, line) == \
                _ReadState.ELECTRONS:
            values = [x.strip() for x in line.split('=')]
            logging.debug('line in electrons section:' + str(line))
            if "mixing_mode" in line:
                mixing_mode = values[1]
                self.mixing_mode = mixing_mode
            elif "mixing_beta" in line:
                mixing_beta = float(values[1])
                self.mixing_beta = mixing_beta
            elif "conv_thr" in line:
                convergence_threshold = values[1]
                # numbers like 1.0d-7 might have to be converted to float
                self.convergence_threshold = convergence_threshold
            line = f.next()
        return line

    def process_ions(self, f):
        print('processing ions section')
        f.next()
        line = f.next()
        return line

    def process_cell(self, f):
        print('processing cell section')
        f.next()
        line = f.next()
        return line

    def process_atomic_species(self, f):
        print('processing atomic species section')
        line = f.next()
        self.SP[CUBA.CHEMICAL_SPECIE] = []
        self.SP[CUBA.MASS] = []
        self.pseudopotential_files = []
        while _ReadState.get_state(_ReadState.ATOMIC_SPECIES, line) == \
                _ReadState.ATOMIC_SPECIES:
            values = line.split()
            logging.debug('line in atomic species section:' + str(line))
#            print('atomtypes:'+str(self.atomtypes))

            if len(values) > 0:
                if values[0] in atomtypes:
                    print("atom type:" + values[0])
                    # self.dc(CHEMICAL_SPECIE = values[0])
                    self.SP[CUBA.CHEMICAL_SPECIE].append(values[0])
                    mass = float(values[1])
                    self.SP[CUBA.MASS].append(mass)
                    potential_file = values[2]
            line = f.next()
        return line

    def process_k_points(self, f, mode='automatic'):
        # skip line
        print('processing k_points section')
        line = f.next()
        self.k_points_mode = mode

        while _ReadState.get_state(_ReadState.K_POINTS, line) == \
                _ReadState.K_POINTS:
            #        print('line:'+str(line))
            values = line.split()
            if len(values):
                K_points = (values)
                print('k points:' + str(K_points))
                self.k_point_values = K_points
            line = f.next()

        return line

    def process_atomic_positions(self, f, units='(angstrom)'):
        print('processing atomic_positions section')
        try:
            line = f.next()
            self.position_units = units
            particle_list = []
            while _ReadState.get_state(_ReadState.ATOMIC_POSITIONS, line) \
                    == _ReadState.ATOMIC_POSITIONS:
                logging.debug('line in atomic positions section:' + str(line))
                values = line.split()
                print('values:' + str(values))
                atom_pos = [0, 0, 0]
                if values[0] in atomtypes:
                    atomtype = values[0]
                    # store position in meters; original in Angstrom
                    atom_pos[0] = float(values[1]) * 1e-10
                    atom_pos[1] = float(values[2]) * 1e-10
                    atom_pos[2] = float(values[3]) * 1e-10
                    # s = str(i)
                    # make uid using string of index
                    # uidstring =
                    # uid = uuid.UUID(uidstring)

                    p = Particle([atom_pos[0], atom_pos[1], atom_pos[2]])
                    # ,uuid.UUID(int=i)
                    print('uid:' + str(p.uid))
                    p.data[CUBA.CHEMICAL_SPECIE] = atomtype
                    particle_list.append(p)
                    try:
                        line = f.next()
                    except StopIteration:
                        print('EOF')
                        break
            if len(particle_list):
                self.pc.add_particles(particle_list)
                #  part_container.add_particle(p)
                n = self.count_particles()
                print('n_particles=' + str(n))
        except StopIteration:
            return
        return

    def count_particles(self,pc=None):
        '''
        Take a particular pc to count, otherwise count all datasets
        '''
        n = 0
        if pc is None:
            names = self.get_dataset_names()
            for name in names:
                n += self.count_particles(self.get_dataset(name))
            return n
        else:
            for particle in pc.iter_particles():
                n += 1
        return n

    def get_particle_masses(self,pc=None):
        '''
        Take a particular pc to count, otherwise count all datasets
        '''
        #todo - take care of possible isotopes - same atom, different mass

        particle_dict = {}
        if pc is None:
            names = self.get_dataset_names()
            for name in names:
                pc = self.get_dataset(name)
                for particle in pc.iter_particles():
                    atomtype_list =particle.data[CUBA.CHEMICAL_SPECIE]
                    atomtype = atomtype_list[0]
                    mass =particle.data[CUBA.MASS]
                    if not atomtype in particle_dict:
                        particle_dict[atomtype]=mass
                        logging.debug('masses:'+str(particle_dict))

            logging.debug('1.atomtypes:'+str(atomtypes))
            n_atom_types = len(atomtypes)
        else:
            #implement case where particular pc is sent
            pass
        return particle_dict

    def count_atom_types(self,pc=None):
        '''
        Take a particular pc to count, otherwise count all datasets
        '''
        atomtypes = set([])
        n_atom_types = 0
        if pc is None:
            names = self.get_dataset_names()
            for name in names:
                n_atom_types += self.count_atom_types(self.get_dataset(name))
            return n_atom_types
        else:
            for particle in pc.iter_particles():
                atomtype_list =particle.data[CUBA.CHEMICAL_SPECIE]
                atomtype = atomtype_list[0]
                logging.debug('0.atomtype:'+str(atomtype))
                atomtypes.add(atomtype)
            logging.debug('1.atomtypes:'+str(atomtypes))
            n_atom_types = len(atomtypes)
        return n_atom_types

    def what_atom_types(self,pc=None):
        atomtype_set = self.what_atom_types_set(pc=pc)
        logging.debug('2.atomtypes:'+str(atomtype_set))
        return [atomtype for atomtype in atomtype_set]


    def what_atom_types_set(self,pc=None):
        '''
        Take a particular pc to count, otherwise count all datasets
        '''
        atomtypes = set([])
        if pc is None:
            names = self.get_dataset_names()
            for name in names:
                more_atoms = self.what_atom_types_set(self.get_dataset(name))
                logging.debug('2.5.more atoms:'+str(more_atoms))
                for atomtype in more_atoms:
                    atomtypes.add(atomtype)
            return atomtypes
        else:
            for particle in pc.iter_particles():
                atomtype_list =particle.data[CUBA.CHEMICAL_SPECIE]
                atomtype = atomtype_list[0]
                logging.debug('3.atomtype:'+str(atomtype))
                atomtypes.add(atomtype)
            logging.debug('4.atomtypes:'+str(atomtypes))
        return atomtypes

    def add_dataset(self, container):
        """Add a CUDS container
        Parameters
        ----------
        container : {ABCParticles}
            The CUDS container to add to the engine.
        Raises
        ------
        TypeError:
            If the container type is not supported (i.e. ABCLattice, ABCMesh).
        ValueError:
            If there is already a dataset with the given name.
        """
        if not isinstance(container, ABCParticles):
            raise TypeError(
                "The type of the dataset container is not supported")

        if container.name in self.get_dataset_names():
            raise ValueError(
                'Particle container \'{}\' already exists'.format(
                    container.name))
        else:
#            self._data_manager.new_particles(container)
#            uname = uuid.uuid4()
            self.datasets[container.name] = container
#            particles = container
#            self._unames[particles.name] = uname
#            self._names[uname] = particles.name

#            lammps_pc = LammpsParticles(self, uname)
#            self._lpcs[uname] = lammps_pc

#            self._handle_new_particles(uname, particles)
#            return lammps_pc

class _ReadState(Enum):
    UNKNOWN, UNSUPPORTED, CONTROL, SYSTEM, ELECTRONS, IONS, CELL, \
        ATOMIC_SPECIES, K_POINTS, ATOMIC_POSITIONS = range(10)

    @staticmethod
    def get_state(current_state, line):
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

def which(name):
    found = 0
    for path in os.getenv("PATH").split(os.path.pathsep):
        full_path = path + os.sep + name
        if os.path.exists(full_path):
            """
            if os.stat(full_path).st_mode & stat.S_IXUSR:
                found = 1
                print(full_path)
            """
            found = 1
            print(full_path)
    # Return a UNIX-style exit code so it can be checked by calling scripts.
    # Programming shortcut to toggle the value of found: 1 => 0, 0 => 1.
    sys.exit(1 - found)

atomtypes = ["C", "H", "He", "N", "O", "Na", "Mg"]

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None



if __name__ == "__main__":
    wrapper = QeWrapper()
#    filename = 'xyzoutput.txt.bak'
#    filename = '../../examples/input_pw.in'
    filename = 'tests/pw.in'
    print('started parsing qe input file ' + str(filename))
    wrapper.read_espresso_input_file(filename)
    print('done parsing qe input file ' + str(filename))

    print('started writing qe input file ' + str(filename))
    new_inputfilename = 'tests/new_input_pw.in'
    wrapper.write_espresso_input_file(new_inputfilename)
    print('done writing qe input file ' + str(new_inputfilename))

    ppfilename = 'tests/testpp.in'
    print('started writing qe pp input file ' + str(ppfilename))
    wrapper.write_espresso_pp_file()
    print('done writing qe pp input file ' + str(ppfilename))

    filename = 'tests/xyzoutput.txt.bak'
    print('started parsing qe output file ' + str(filename))
    wrapper.read_espresso_output_file(filename)
    print('finished parsing qe output file ' + str(filename))
