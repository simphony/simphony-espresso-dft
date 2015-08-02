import string
from enum import Enum
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer
from simphony.cuds.particles import Particle, Particles


def WriteEspressoInputFile(file_name,data_container,particle_container):
    """
    :param file_name: name of the input file to write
    :param data_container: contains data to write to file
    :return:
    """
 #   d=DataContainer.items()
 #   d2=data_container.items()
    SP = data_container
    pc = particle_container
    print('attempting to write '+file_name)
    with open(file_name, 'w') as f:
        try:
            #CONTROL section
            #apparently a comma is not required at the end of every line
            line = '&CONTROL\n'  #calculation
            f.write(line)
            if CUBA.TORQUE in SP:
                line = '\t calculation=\''+str(SP[CUBA.TORQUE])+'\'\n'   #calculation
                f.write(line)
            if CUBA.ZETA_POTENTIAL in SP:
                line = '\t restart_mode=\''+str(SP[CUBA.ZETA_POTENTIAL])+'\'\n' #restart_mode
                f.write(line)
            if CUBA.YOUNG_MODULUS in SP:
                line = '\t pseudo_dir=\''+str(SP[CUBA.YOUNG_MODULUS])+'\'\n'  #pseudo dir
                f.write(line)
            if CUBA.VOLUME_FRACTION in SP:
                line = '\t prefix=\''+str(SP[CUBA.VOLUME_FRACTION])+'\'\n'  #prefix
                f.write(line)
            if CUBA.AMPHIPHILICITY in SP:
                line = '\t tprnfor='+str(SP[CUBA.AMPHIPHILICITY])+'\n'  #tprnfor
                f.write(line)
            if CUBA.NUMBER_OF_TIME_STEPS in SP:
                line = '\t max_seconds='+str(SP[CUBA.NUMBER_OF_TIME_STEPS])+'\n'
                f.write(line)
            if CUBA.DIRECTION in SP:
                line = '\t outdir=\''+str(SP[CUBA.DIRECTION])+'\'\n'  #outdir
                f.write(line)
            line = '/\n'
            f.write(line)

            #SYSTEM section
            line = '&SYSTEM\n'  #calculation
            f.write(line)
            if CUBA.ROLLING_FRICTION in SP:
                line = '\t ibrav='+str(SP[CUBA.ROLLING_FRICTION])+'\n'  #outdir
                f.write(line)
            if CUBA.ORIGINAL_POSITION in SP:
                line = '\t celldm(1)='+str(SP[CUBA.ORIGINAL_POSITION][0])+'\n'
                f.write(line)
                line = '\t celldm(2)='+str(SP[CUBA.ORIGINAL_POSITION][1])+'\n'
                f.write(line)
                line = '\t celldm(3)='+str(SP[CUBA.ORIGINAL_POSITION][2])+'\n'
                f.write(line)

             # here goes nat and ntype


            if CUBA.LN_OF_RESTITUTION_COEFFICIENT in SP:
                line = '\t ecutwfc='+str(SP[CUBA.LN_OF_RESTITUTION_COEFFICIENT])+'\n'
                f.write(line)
            if CUBA.POISSON_RATIO in SP:
                line = '\t ecutrho='+str(SP[CUBA.POISSON_RATIO])+'\n'
                f.write(line)
            if CUBA.LATTICE_SPACING in SP:
                line = '\t input_dft=\''+str(SP[CUBA.LATTICE_SPACING])+'\'\n'  #outdir
                f.write(line)
            line = '/\n'
            f.write(line)


            #ELECTRONS section
            line = '&ELECTRONS\n'  #calculation
            f.write(line)
            if CUBA.SMOOTHING_LENGTH in SP:
                line = '\t mixing_mode=\''+str(SP[CUBA.SMOOTHING_LENGTH])+'\'\n'
                f.write(line)
            if CUBA.PHASE_INTERACTION_STRENGTH in SP:
                line = '\t mixing_beta='+str(SP[CUBA.PHASE_INTERACTION_STRENGTH])+'\n'
                f.write(line)
            if CUBA.DEBYE_LENGTH in SP:
                line = '\t conv_thr='+str(SP[CUBA.DEBYE_LENGTH])+'\n'
                f.write(line)
            line = '/\n'
            f.write(line)

            #IONS section
            line = '&IONS\n'  #calculation
            f.write(line)
            line = '/\n'
            f.write(line)

            #CELL section
            line = '&CELL\n'  #calculation
            f.write(line)
            line = '/\n'
            f.write(line)

            #CELL section
            line = 'ATOMIC_SPECIES\n'  #calculation
            f.write(line)

            if CUBA.CHEMICAL_SPECIE in SP:
                for i in range(0,len(SP[CUBA.CHEMICAL_SPECIE])):
                    line = '\t'+str(SP[CUBA.CHEMICAL_SPECIE][i])+' '+\
                           str(SP[CUBA.MASS][i])+' \''+str(SP[CUBA.FRICTION_COEFFICIENT][i])+'\'\n'
                    f.write(line)
            line = '\n'
            f.write(line)

            if CUBA.PROBABILITY_COEFFICIENT in SP:
                line ='K_POINTS '+str(SP[CUBA.PROBABILITY_COEFFICIENT])+'\n'
                f.write(line)
            if CUBA.EQUATION_OF_STATE_COEFFICIENT in SP:
                line = ''
                for k_point in SP[CUBA.EQUATION_OF_STATE_COEFFICIENT]:
                    line = line+str(k_point)+' '
                f.write(line)
            f.write('\n\n')

            if CUBA.KINEMATIC_VISCOSITY in SP:
                line ='ATOMIC_POSITIONS '+str(SP[CUBA.KINEMATIC_VISCOSITY])+'\n'
                f.write(line)

            #this is until setting/getting particles from particle_container works
            kluge(f)
            return

            for particle in pc.iter_particles:
                atom_type = particle.data
                atom_coords = [particle.coordinates[0],
                               particle.coordinates[1],
                               particle.coordinates[2]]
                line =str(atom_type)+' '+str(particle.coordinates[0])+' '+\
                      str(particle.coordinates[1])+' '+str(particle.coordinates[2])
                f.write(line)


        except:
                ('error in write block of WriteEspressoInputFile')
                raise
    f.closed

def kluge(f):
    line =' C 1.0 2.0 3.0 \n'
    f.write(line)
    line =' C 3.0 3.0 4.0 \n'
    f.write(line)
    line =' C 3.0 4.0 5.0 \n'
    f.write(line)
    line =' C 4.0 5.0 6.0 \n'
    f.write(line)
    return




if __name__ == "__main__":
    filename = '../../examples/input_pw.in'
    print('started parsing file '+str(filename))
    ReadEspressoInputFile(filename)
