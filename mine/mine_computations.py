from mine_utils import *
import numpy as np
from math import log

def entropy(probs): 
    return -sum(p * log(p, 2) for p in probs if p > 0)
                 
def I(P_ordinals, Q_map):
    return HQ(Q_map) + HP(P_ordinals) - HPQ(P_ordinals, Q_map)

def GetProbabilityDistribution(P):
    partittions = GroupPointsByPartition(P)
    # The probability mass of each partition is the fraction of points that lie in this partition
    prob_mass = lambda p: len(partittions[p]) / float(len(P))
    return map(prob_mass, partittions)
def test_H():
    assert H(P=[0.25, 0.25, 0.25, 0.25]) == 2
    
    # OpenMIC test case
    assert H(P=[1. / 8, 1. / 4, 1. / 8, 1. / 2]) == 7. / 4
    
    # Test case for joint partition
    """
      4  |   |     |x
      3  |   |  x  |
      2  |   |x    |
         *----+---+-----+-
      1  |x x|     |
         *----+---+-----+-
      0 x|   |    x|
        0|1 2|3 4 5 6
         |   |     |
     """
    # Constructing the grid above
    P = {(0, 0): 0,
         (1, 1): 1,
         (2, 1): 1,
         (3, 2): 2,
         (4, 3): 2,
         (5, 0): 2,
         (6, 4): 3
         }
    
    Q = {(0, 0): 0,
         (1, 1): 1,
         (2, 1): 1,
         (3, 2): 2,
         (4, 3): 2,
         (5, 0): 0,
         (6, 4): 2
         }
    
    # visualize(P, Q)
    
    # Joint entropy computation
    assert (H(P, Q) == H(
                         [0  , 0    , 2. / 7 , 1. / 7,
                          0  , 2. / 7 , 0. / 7 , 0   ,
                        1. / 7 , 0    , 1. / 7 , 0   ]))
     
def HP(Dx, P_ordinals):
    assert is_sorted_increasing_by(Dx, 'x')
    
    #Number of points in the partition
    m = P_ordinals[-1] + abs(P_ordinals[0])
    return entropy(np.array(GetPartitionHistogram(Dx, P_ordinals)) / float(m))

def HQ(Q_map):
    pass

def HPQ(P_ordinals, Q_map):
    pass
