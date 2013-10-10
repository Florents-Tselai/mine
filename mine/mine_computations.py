from mine_utils import *
import numpy as np
from math import log

def entropy(probs): 
    return -sum(p * log(p, 2) for p in probs if p > 0)
                 
def I(P_ordinals, Q_map):
    return HQ(Q_map) + HP(P_ordinals) - HPQ(P_ordinals, Q_map)

def HP(Dx, P_ordinals):
    assert is_sorted_increasing_by(Dx, 'x')
    
    #Number of points in the partition
    m = P_ordinals[-1] + abs(P_ordinals[0])
    return entropy(np.array(GetPartitionHistogram(Dx, P_ordinals)) / float(m))

def HQ(Q_map):
    n = len(Q_map)
    temp = GroupPointsByPartition(Q_map)
    return entropy([len(temp[k]) / float(n) for k in temp])

def HPQ(P_ordinals, Q_map):
    Dx = sort_D_increasing_by(Q_map.keys(), 'x')
    
    x_axis_partition = GetPartitionMapFromOrdinals(Dx, P_ordinals, axis='x')
    m = len(x_axis_partition.keys())
    
    return entropy(np.array(GetGridHistogram(Q_map, x_axis_partition)) / float(m))