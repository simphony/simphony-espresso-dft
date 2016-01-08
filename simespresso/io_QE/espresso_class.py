import logging
import subprocess
import sys

from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_particles import ABCParticles

import data_handler

logging.basicConfig(level=logging.DEBUG)


class QeFileIO(ABCModelingEngine):
    '''
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
    pwname = self.datahandler.input_pwname
        self.datahandler.write_espresso_input_file(pwname)
    path = self.datahandler.path_to_espresso
        logging.debug('path to espresso:'+path)
        if not which(self.path_to_espresso):
            logging.debug('no path to espresso')
            raise ValueError(
                'espresso command not found (looking for '
                            + self.path_to_espresso+')')
        if self.datahandler.mpi:
            command = 'mpirun -np ' + str(self.mpi_Nprocessors) + ' ' + \
                     path + ' < ' + pwname + ' > ' \
                      + self.datahandler.output_filename
        else:
            command = self.path_to_espresso + ' < ' + pwname + ' > ' \
                      + self.datahandler.output_filename
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
            self.datahandler.datasets[container.name] = container

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
            del self.datahandler.datasets[dataset.name]

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
            dataset = self.datahandler.datasets[dataset_name]
            return dataset

    def get_dataset_names(self):
        """Get names of datasets in container
        """
        dataset_names = [name for name in self.datahandler.datasets]
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

