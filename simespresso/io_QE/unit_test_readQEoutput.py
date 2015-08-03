__author__ = 'jeremy'

import unittest
import espresso_output_file_read


class OutcomesTest(unittest.TestCase):

    def test_espresso_data_file_read(self):
        filename = 'xyzoutput.txt'
        print('started parsing file '+str(filename))
        espresso_output_file_read.ReadEspressoOutputFile(filename)

if __name__ == '__main__':
    unittest.main()

