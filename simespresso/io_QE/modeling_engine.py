from abc import ABCMeta, abstractmethod
from espresso_class import qe_functions

class ABCModelingEngineForQE(object):  # pragma: no cover
    """Abstract base class for modeling engines in SimPhoNy.
    Through this interface, the user controls and interacts with the
    simulation/calculation (which is being performed by the modeling
    engine).
    Attributes
    ----------
    BC : DataContainer
        container of attributes related to the boundary conditions
    CM : DataContainer
        container of attributes related to the computational method
    SP : DataContainer
        container of attributes related to the system parameters/conditions
    materials : Materials
        materials related to state data
    """

    def run(self):
        """ Run the modeling engine
        Run the modeling engine using the configured settings (e.g. CM, BC,
        and SP) and the configured state data (e.g. particle, mesh and
        lattice data).
        """
        name_in = './pwtest.in'
        name_out = './pwtest.out'
        mpi = False
        mpi_Nprocessors = 2

        qe_functions.start_qe(self,name_in,name_out,mpi=mpi,mpi_Nprocessors=mpi_Nprocessors)

    def add_dataset(self, container):
        """Add a CUDS container
        Parameters
        ----------
        container : {ABCMesh, ABCParticles, ABCLattice}
            The CUDS container to add to the engine.
        Raises
        ------
        TypeError:
            If the container type is not supported by the engine.
        ValueError:
            If there is already a dataset with the given name.
        """

    def remove_dataset(self, name):
        """ Remove a dataset from the engine
        Parameters
        ----------
        name: str
            name of CUDS container to be deleted
        Raises
        ------
        ValueError:
            If there is no dataset with the given name
        """

    def get_dataset(self, name):
        """ Get the dataset
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

    def get_dataset_names(self):  # pragma: no cover
        """ Returns a list of the datasets' names in the engine workspace.
        """

    def iter_datasets(self, names=None):  # pragma: no cover
        """ Returns an iterator over a subset or all of the containers.
        Parameters
        ----------
        names : sequence of str, optional
            names of specific containers to be iterated over. If names is not
            given, then all containers will be iterated over.
        """