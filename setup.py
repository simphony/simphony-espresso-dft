import os

from setuptools import setup, find_packages

with open('README.rst', 'r') as readme:
    README_TEXT = readme.read()

VERSION = '0.0.1.dev0'


def write_version_py(filename=None):
    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__), 'simespresso', 'version.py')
    ver = """\
version = '%s'
"""
    fh = open(filename, 'wb')
    try:
        fh.write(ver % VERSION)
    finally:
        fh.close()


write_version_py()

setup(
    name='simespresso',
    version=VERSION,
    author='SimPhoNy, EU FP7 Project (Nr. 604005) www.simphony-project.eu',
    description='The espresso wrapper for the SimPhoNy framework',
    long_description=README_TEXT,
    entry_points={
        'simphony.engine': ['espresso = simespresso']},
    packages=find_packages(),
    )
