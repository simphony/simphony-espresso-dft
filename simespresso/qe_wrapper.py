""" qe SimPhoNy Wrapper
This module provides a wrapper for  quantum-espresso
"""
import contextlib
import os
import shutil
import tempfile

from simphony.core.data_container import DataContainer
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_particles import ABCParticles

from simespresso.io_QE.espresso_fileio_data_manager import QeFileIoDataManager
from simespresso.io_QE.qe_process import QeProcess
from qeCubaExtensions import qeCUBAExtension

import logging

@contextlib.contextmanager
def _temp_directory():
    """ context manager that provides temp directory
    The name of the created temp directory is returned when context is entered
    and this directory is deleted when context is exited
    """
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    #shutil.rmtree(temp_dir)


class QeWrapper(ABCModelingEngine):
    """ Wrapper to quantum espresso
    """
    def __init__(self):
        """ Constructor.
        Parameters
        ----------
        """
        self._executable_name = "/home/jr/sw/espresso-5.2.1/bin/pw.x"
        self._executable_name = "pw.x"
        self.BC = DataContainer()
        self.CM = DataContainer()
        self.SP = DataContainer()
        self.SD = DataContainer()
        self.CM_extension = {}  #defaults below

        self.SP_extension = {}
        #self.SP_extension[qeCUBAExtension.PSEUDO_POTENTIAL] = 'vdw-df-c09' #default

        self.BC_extension = {}
        self._data_manager = QeFileIoDataManager(self)


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
        # todo-what container types do we need to support for QE
        if not isinstance(container, ABCParticles):
            raise TypeError(
                "The type of the dataset container is not supported")

        logging.debug('adding container '+str(container.name))
        if container.name in self._data_manager._pc_cache:
            raise ValueError(
                'Particle container \'{}\' already exists'.format(
                    container.name))
        else:
            self._data_manager.add_particles(container)

    def get_dataset(self, name):
        """ Get the dataset
        The returned particle container can be used to query
        and change the related data inside LAMMPS.
        Parameters
        ----------
        name: str
            name of CUDS container to be retrieved.
        Returns
        -------
        container :
            A proxy of the dataset named ``name`` that is stored
            internally in the Engine.
        Raises
        ------
        ValueError:
            If there is no dataset with the given name
        """
        if name in self._data_manager._pc_cache:
            return self._data_manager._pc_cache[name]
        else:
            raise ValueError(
                'Particle container \'{}\` does not exist'.format(name))

    def get_dataset_names(self):
        """ Returns the names of all the datasets
        """
        # TODO  (simphony-common #218)
        return [name for name in self._data_manager._pc_cache]

    def remove_dataset(self, name):
        """ Remove a dataset
        Parameters
        ----------
        name: str
            name of CUDS container to be deleted
        Raises
        ------
        ValueError:
            If there is no dataset with the given name
        """
        if name in self._data_manager._pc_cache:
            del self._data_manager._pc_cache[name]
        else:
            raise ValueError(
                'Particles \'{}\' does not exist'.format(name))

    def iter_datasets(self, names=None):
        """ Returns an iterator over a subset or all of the containers.
        Parameters
        ----------
        names : sequence of str, optional
            names of specific containers to be iterated over. If names is not
            given, then all containers will be iterated over.
        """
        if names is None:
            for name in self._data_manager:
                yield self._data_manager[name]
        else:
            for name in names:
                if name in self._data_manager:
                    yield self._data_manager[name]
                else:
                    raise ValueError(
                        'Particle container \'{}\` does not exist'.format(
                            name))

    def run(self):
        """ Run qe based on configuration and data
        """
        with _temp_directory() as temp_dir:
            input_data_filename = os.path.join(
                temp_dir, "data_in.qe")
            output_data_filename = os.path.join(
                temp_dir, "data_out.qe")

            # before running, we flush any changes
            self._data_manager.flush(input_data_filename)

#            commands = self._script_writer.get_configuration(
 #               input_data_file=input_data_filename,
  #              output_data_file=output_data_filename,
   #             BC=_combine(self.BC, self.BC_extension),
    #            CM=_combine(self.CM, self.CM_extension),
     #           SP=_combine(self.SP, self.SP_extension))

            BC=_combine(self.BC, self.BC_extension)
            CM=_combine(self.CM, self.CM_extension)
            SP=_combine(self.SP, self.SP_extension)


            process = QeProcess(self._data_manager, qe_executable=self._executable_name,
                                    log_directory=temp_dir)
            process.run(input_data_filename,output_data_filename,BC,CM,SP)

            # after running, we read any changes from lammps
            self._data_manager._read_espresso_output_file(output_data_filename)
    #        self._data_manager._read_densities(n_latticepoints, file_iter, aviz_filename=False):



def _combine(data_container, data_container_extension):
    """ Combine a the approved CUBA with non-approved CUBA key-values
    Parameters
    ----------
    data_container : DataContainer
        data container with CUBA attributes
    data_container_extension : dict
        data container with non-approved CUBA attributes
    Returns
    ----------
    dict
        dictionary containing the approved adn non-approved
        CUBA key-values
    """
    result = dict(data_container_extension)
    result.update(data_container) #combines result,data_container to 1 dict
    return result

