from enum import Enum, unique

@unique
class CUBAExtension(Enum):
    """ Provisional CUBA keywords specific for Quantum Espresso

    These are additional CUBA-Keywords that are not included
    in simphony-common yet. The proposed description for
    these new CUBA keywords is:

    - description: Simulation box faces
    domain: [MD]
    key: BOX_FACES
    name: BoxFaces
    number: 100
    shape: [1]
    type: double
    - description: Simulation box origin
    domain: [MD]
    key: BOX_ORIGIN
    name: BoxOrigin
    number: 102
    shape: [3]
    type: double
    - description: Espresso simulation type
    domain: [MD]
    key: SIMULATION_TYPE
    name: SimulationType
    number: 102
    shape: [1]
    type: string

"""

    BOX_FACES = "BOX_FACES"
    BOX_VECTORS = "BOX_VECTORS"
    BOX_ORIGIN = "BOX_ORIGIN"
    SIMULATION_TYPE = "SIMULATION_TYPE"