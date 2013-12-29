import unittest
from computations import *
import numpy as np
from numpy import array
from numpy.testing import assert_array_equal

class computations_test(unittest.TestCase):
    def setUp(self):        
        #Example from wiki https://en.wikipedia.org/wiki/Marginal_distribution
        # Appears also in Cover et.al. book pg. 18
        self.wiki_example = array([[4./32, 2./32, 1./32, 1./32], 
                                   [2./32, 4./32, 1./32, 1./32], 
                                   [2./32, 2./32, 2./32, 2./32], 
                                   [8./32, 0., 0., 0.]])
        
        # Example from Guttler's presentation
        self.guttler_example = array([[1./20, 1./20, 1./20], 
                                      [3./20, 5./20, 2./20], 
                                      [3./20, 2./20, 2./20]])

        assert np.sum(self.wiki_example) == 1.
        
    def test_H(self):
        X = getXDistribution(self.wiki_example)
        assert H(X) == 7./4
        Y = getYDistribution(self.wiki_example)
        assert H(Y) == 2.
        
        #Joint entropy
        assert H(self.wiki_example) == 27./8
        
    def test_getXDistribution(self):
        assert_array_equal(getXDistribution(self.wiki_example), array([1./2, 1./4, 1./8, 1./8]))
        
    def test_getYDistribution(self):
        assert_array_equal(getYDistribution(self.wiki_example), array([1./4, 1./4, 1./4, 1./4]))
        
    def test_I(self):
        assert I(self.wiki_example) == 3./8
        

if __name__ == '__main__':
    unittest.main()