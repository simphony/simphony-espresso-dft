import unittest

from simespresso.io_QE import espresso_class


class OutcomesTest(unittest.TestCase):

    def test_start_qe(self):
        print('testing the starting of qe')
        name_in = './pwtest.in'
        name_out = './pwtest.out'
        mpi = False
        mpi_Nprocessors = 2
        path_to_espresso = 'pw.x'
        wrapper = espresso_class.qe_functions()
        wrapper.start_qe(name_in, name_out,
                         path_to_espresso=path_to_espresso,
                         mpi=mpi, mpi_Nprocessors=mpi_Nprocessors)

if __name__ == '__main__':
    unittest.main()
