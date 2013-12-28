import unittest
from computations import *
import numpy as np
from numpy import array
from numpy.testing import assert_array_equal

class computations_test(unittest.TestCase):
    def setUp(self):        
        #Example from wiki https://en.wikipedia.org/wiki/Marginal_distribution
        # Appears also in Cover et.al. book pg. 18
        self.example = array([[4./32, 2./32, 1./32, 1./32], 
                                 [2./32, 4./32, 1./32, 1./32], 
                                 [2./32, 2./32, 2./32, 2./32], 
                                 [8./32, 0., 0., 0.]])

        assert np.sum(self.example) == 1.
        
    def test_H(self):
        X = getXDistribution(self.example)
        assert H(X) == 7./4
        Y = getYDistribution(self.example)
        assert H(Y) == 2.
        
        #Joint distribution
        assert H(self.example) == 27./8
        
    def test_getXDistribution(self):
        X = np.sum(self.example, axis=0)
        assert_array_equal(X, array([1./2, 1./4, 1./8, 1./8]))
        
    def test_getYDistribution(self):
        Y = np.sum(self.example, axis=1)
        assert_array_equal(Y, array([1./4, 1./4, 1./4, 1./4]))

if __name__ == '__main__':
    unittest.main()