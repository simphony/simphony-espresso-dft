from enum import Enum
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.particles import Particle, Particles
from simphony.cuds.lattice import Lattice
import numpy as np
#from simphony.cuds.abstractparticles import ABCParticles

import uuid
import logging
logging.basicConfig(level=logging.DEBUG)



class qe_wrapper(object):
    '''
    functions for reading and writing quantum espresso input and output files
    '''

    def __init__(self):
        self.SP = DataContainer()   #System Parameters and Conditions


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

                    L = Lattice('quantum espresso lattice',Ltype,celldm,n_lattice_points,[0,0,0])

                    #4th line - don't care
                    line = file_iter.next()
                    logging.debug('read line4:'+str(line))
                    values = line.split()
                    floats = [float(val) for val in values]

            except StopIteration:
                print('eof reached')
                logging.warning('eof reached')
            except :
                print("problem with line", line)
                logging.warning("problem with line", line)
                raise
            #TODO skip all the atom definition lines, eg look for alphabetic characters at beginning of line
            logging.debug('skipping '+str(n_atoms)+' lines')
            for i in range(0,n_atoms+1):
                    line = file_iter.next()

            logging.debug('read '+str(n_lattice_points)+' lattice point lines')
            L = self.read_densities(n_lattice_points,L,file_iter,aviz_filename='avizout.xyz')
            #put pc into dc!



    def running_index_to_node_index(self,index,n_latticepoints):
        node_z = index/(n_latticepoints[0]*n_latticepoints[1])
        node_y = (index/n_latticepoints[0])%n_latticepoints[1]
        node_x = index%n_latticepoints[0]
        return [node_x,node_y,node_z]

    def read_densities(self,n_latticepoints,L,file_iter,aviz_filename=False):
        charge_density = self.read_xyz(n_latticepoints,file_iter)
        charge_as_list = charge_density.flatten()
        nodelist = []
        for charge_index, charge in enumerate(charge_as_list):
            node_index = self.running_index_to_node_index(charge_index,n_latticepoints)# TODO convert charge_index to node index
            node = L.get_node(node_index)
            node.data[CUBA.MASS] = charge
            nodelist.append(node)
        L.update_nodes(nodelist)  #we have to ask the lattice to update the changed node


        if aviz_filename:
            self.write_aviz_output(charge_density,aviz_filename,L.base_vect)
        return L

    #    iterator = L.iter_nodes(n_latticepoints)
    #    for i in range(0,n_latticepoints[0]):
    #        for j in range(0,n_latticepoints[1]):
    #            for k in range(0,n_latticepoints[2]):
    #                L.data = charge_density[i,j,k]
     #               iterator.next()
    def write_aviz_output(self,xyz_array,aviz_xyzfile,base_vector):
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
            logging.debug('Got fewer points than expected')
        return charge_density

if __name__ == "__main__":
    filename = 'xyzoutput.txt.bak'
    filename = 'xyzoutput.txt'
    print('started parsing file '+str(filename))
    wrapper = qe_wrapper()
    wrapper.ReadEspressoOutputFile(filename)
