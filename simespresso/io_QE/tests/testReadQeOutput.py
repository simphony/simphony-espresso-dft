__author__ = 'jeremy'


import unittest
import os
from simespresso.io_QE import espresso_class
import logging

class OutcomesTest(unittest.TestCase):

    def test_read_espresso_output_file(self):
        file_name = 'pwtest.out'
        if not(os.path.exists(file_name)):
            import logging
         #   logging.debug("file "+str(file_name)+" not found")
            print("file "+str(file_name)+" not found")
            return(1)
        print('testing reading of qe output file '+str(file_name))
        qe_wrapper = espresso_class.qe_functions()
        qe_wrapper.read_espresso_output_file(file_name)


    def test_running_index_to_node_inredex(self):
        print('testing espresso_data_file_read')
        n_latticepoints = [10,7,6]
        index = 5
        espresso_wrapper = espresso_class.qe_functions()
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

