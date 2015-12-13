from enum import Enum, unique


@unique
class CUBAExtension(Enum):
    """ Provisional CUBA keywords specific for Lammps-Md
    These are additional CUBA-Keywords that are not included
    in simphony-common yet. The proposed description for
    these new CUBA keywords is:
    - description: Simulation box vectors
    domain: [MD]
    key: BOX_VECTORS
    name: BoxVectors
    number: 101
    shape: [3,3]
    type: double
    - description: Simulation box origin
    domain: [MD]
    key: BOX_ORIGIN
    name: BoxOrigin
    number: 102
    shape: [3]
    type: double
    - description: Pair potentials
    domain: [MD]
    key: PAIR_POTENTIALS
    name: PairPotentials
    number: 104
    shape: [20]
    type: string
"""

    BOX_VECTORS = "BOX_VECTORS"
    BOX_ORIGIN = "BOX_ORIGIN"
    PAIR_POTENTIALS = "PAIR_POTENTIALS"