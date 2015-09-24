__author__ = 'jeremy'

import unittest
import espresso_wrapper


class OutcomesTest(unittest.TestCase):

    def test_ReadEspressoOutputFile(self):
        filename = 'xyzoutput.txt'
        print('testing parsing file '+str(filename))
        espresso_wrapper.ReadEspressoOutputFile(filename)


    def test_running_index_to_node_index(self):
        print('testing espresso_data_file_read')
        n_latticepoints = [10,7,6]
        index = 5
        indices = espresso_wrapper.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[5,0,0])

        index = 15
        indices = espresso_wrapper.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[5,1,0])

        index = 23
        indices = espresso_wrapper.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[3,2,0])

        index = 70
        indices = espresso_wrapper.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[0,0,1])

        index = 75
        indices = espresso_wrapper.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[5,0,1])

        index = 150
        indices = espresso_wrapper.running_index_to_node_index(index,n_latticepoints)
        print('index '+str(index)+' indices:'+str(indices))
        self.assertTrue(indices==[0,1,2])



if __name__ == '__main__':
    unittest.main()

