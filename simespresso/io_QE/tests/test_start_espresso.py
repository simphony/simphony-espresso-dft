
import unittest
import os
from simespresso.io_QE import espresso_class


class OutcomesTest(unittest.TestCase):

    def test_start_qe(self):
        name_in = './test_pw.in'
        name_out = './test_pw.out'
        mpi=False
        mpi_Nprocessors=2
        path_to_espresso = '/usr/bin/pw.x'
        print('qe wrapper attempting to run: '+command)
        start_qe(self,name_in,name_out,path_to_espresso=path_to_espresso,mpi=mpi,mpi_Nprocessors=mpi_Nprocessors):


if __name__ == '__main__':
    unittest.main()



