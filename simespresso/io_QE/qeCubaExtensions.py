from enum import Enum, unique


@unique
class qeCUBAExtension(Enum):
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
    - description: Simulation box faces
    domain: [MD]
    key: BOX_FACES
    name: BoxFaces
    number: 100
    shape: [1]
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
    - description: Pseudopotentials
    domain: [MD]
    key: PAIR_POTENTIALS
    name: PairPotentials
    number: 200
    shape: [20]
    type: string
    - description: kpoint sampling method
    domain: [MD]
    key: K_POINT_SAMPLING_METHOD
    name: K_point_sampling_method
    number: 210
    shape: [20]
    type: string
    - description: k point sampling
    domain: [MD]
    key: K_POINT_SAMPLING
    name: K_point_sampling
    number:220
    shape: [3]
    type: double
    - description: total energy
    domain: [MD]
    key: TOTAL_ENERGY
    name: TOTAL_ENERGY
    number:1234
    shape: 1
    type: double
    """

    BOX_VECTORS = "BOX_VECTORS"
    BOX_FACES = "BOX_FACES"
    BOX_ORIGIN = "BOX_ORIGIN"
    PAIR_POTENTIALS = "PAIR_POTENTIALS"
    PSEUDO_POTENTIAL = "PSEUDO_POTENTIAL"
    K_POINT_SAMPLING_METHOD = "K_POINT_SAMPLING_METHOD"
    K_POINT_SAMPLING = "K_POINT_SAMPLING"
    TOTAL_ENERGY = "TOTAL_ENERGY"

    # add 'extra results field with total energy, iter, delta E
