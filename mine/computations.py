import numpy as np
from sklearn.metrics import mutual_info_score

EPSILON = 1e-8

def H(distribution):
    def entropy(P):
        '''
        Return the Shannon entropy of a probability vector P
        See http://www.scholarpedia.org/article/Entropy#Shannon_entropy
        '''
        h = -sum(i*np.log2(i) for i in P if i>0)
        assert h >= 0
        assert h <= np.log2(len(P)) + EPSILON
        return h
    # Flatter a 2-dim (grid histogram array)
    return entropy(distribution.ravel())

def getXDistribution(grid_histogram):
    return np.sum(grid_histogram, axis=0)

def getYDistribution(grid_histogram):
    return np.sum(grid_histogram, axis=1)

def I(joint_distribution_histogram):
    x_distribution = getXDistribution(joint_distribution_histogram)
    y_distribution = getYDistribution(joint_distribution_histogram)
    joint_distribution = joint_distribution_histogram.ravel()
    return H(x_distribution) + H(y_distribution) - H(joint_distribution)

