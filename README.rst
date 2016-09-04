simphony-espresso
===============

The quantum espresso engine-wrapper for the SimPhoNy framework (www.simphony-project.eu).

.. image:: https://api.travis-ci.org/simphony/simphony-espresso-dft.svg?branch=master
   :target: https://travis-ci.org/simphony/simphony-espresso-dft
   :alt: Build status

.. image:: https://codecov.io/github/simphony/simphony-espresso-dft/coverage.svg?branch=io6
   :target: http://codecov.io/github/simphony/simphony-espresso-dft/?branch=master
   :alt: Test coverage

.. image:: https://readthedocs.org/projects/simphony-espresso-dft/badge/?version=master
   :target: https://readthedocs.org/projects/simphony-espresso-dft/?badge=master
   :alt: Documentation Status


Repository
----------

simphony-espresso is hosted on github: https://github.com/simphony/simphony-espresso-dft

Requirements
------------

- pyyaml >= 3.11
- `simphony-common`_ >= 0.2.0

Optional requirements
~~~~~~~~~~~~~~~~~~~~~

TTo take advantage of mpi it needs to be installed, for example

- sudo apt-get install libcr-dev mpich2 mpich2-doc

To support the documentation built you need the following packages:

- sphinx >= 1.2.3
- sphinxcontrib-napoleon >= 0.2.10

Installation
------------

The package requires python 2.7.x. Installation is based on setuptools::

    # build and install
    python setup.py install

or::

    # build for in-place development
    python setup.py develop

quantum espresso installation
~~~~~~~~~~~~~~~~~~~

This engine-wrapper uses the quantum espresso simulation engine, which on ubuntu linux may be installed with

    sudo apt-get install build-essential fftw3-dev gfortran
    sudo apt-get install quantum-espresso

For general quantum-espresso install information, see http://www.quantum-espresso.org/wp-content/uploads/Doc/user_guide.pdf


Usage
-----

After installation, the user should be able to import the ``liggghts`` engine plugin module by::

  from simphony.engine import simespresso

    engine = simespresso.qe_wrapper()


Testing
-------

To run the full test-suite run::

    python -m unittest discover

Directory structure
-------------------

- simespresso -- holds the simespresso wrapper implementation
    - io -- file-io related communication with LIGGGHTS
        - tests -- testing related utilities
- examples -- holds different examples

.. _simphony-common: https://github.com/simphony/simphony-common