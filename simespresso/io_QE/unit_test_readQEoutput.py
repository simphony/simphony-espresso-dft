__author__ = 'jeremy'

import unittest
import espresso_output_file_read


class OutcomesTest(unittest.TestCase):

    def test_espresso_data_file_read(self):
        filename = 'xyzoutput.txt'
        print('started parsing file '+str(filename))
        espresso_output_file_read.ReadEspressoOutputFile(filename)


    def test_espresso_data_file_read(self):
        n_latticepoints = [10,7,6]
        index = 5
        indices = espresso_output_file_read.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[5,0,0])

        index = 15
        indices = espresso_output_file_read.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[5,1,0])

        index = 23
        indices = espresso_output_file_read.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[3,2,0])

        index = 70
        indices = espresso_output_file_read.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[0,0,1])

        index = 75
        indices = espresso_output_file_read.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[5,0,1])

        index = 150
        indices = espresso_output_file_read.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[0,1,2])



if __name__ == '__main__':
    unittest.main()

