from collections import OrderedDict
from mine import visualize_partition, EquipartitionYAxis, GetClumpsPartition, H, GetPartitionIndices, GroupPartitionsPoints, GetProbabilityDistribution
import numpy as np

"""
Example from Albanese et. al. pg. 4
http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf

"""

D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
x_partition = [1.8, 2.2, 7.8, 8.2]
y_partition = [2.5, 4.6]

# visualize_partition(D, x_partition, y_partition)

def test_EquipartitionYAxis():
    #Albanese data
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
    D = sorted(D, key=lambda p: p[1])
    Q = EquipartitionYAxis(D, y=3)
    
    assert Q[(1, 1)] == 1
    assert Q[(1, 2)] == 1
    assert Q[(1, 3)] == 2
    assert Q[(1, 4)] == 2
    assert Q[(2, 3)] == 2
    assert Q[(2, 4)] == 2
    assert Q[(3, 5)] == 3
    assert Q[(4, 6)] == 3
    assert Q[(5, 6)] == 3
    assert Q[(6, 6)] == 3
    assert Q[(7, 5)] == 3
    assert Q[(8, 3)] == 2
    assert Q[(9, 2)] == 1
    assert Q[(9, 1)] == 1
    
    #Spinellis OpenMIC data
    D = [(0, 0), (1, 1), (3, 2), (2, 1), (5, 0), (4, 3), (6, 4)]
    Q = EquipartitionYAxis(D, y=3)
    partition_groups = GroupPartitionsPoints(Q)
    assert (0, 0) and (5, 0) in partition_groups[1]
    assert (1, 1) and (2, 1) in partition_groups[2]
    assert (3, 2) and (4, 3) and (6, 4) in partition_groups[3]
    
def test_GetClumpsPartition():
    
    # Albanese et.al. data
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]

    Q = {}
    Q[(1, 1)] = 1
    Q[(1, 2)] = 1
    Q[(1, 3)] = 2
    Q[(1, 4)] = 2
    Q[(2, 3)] = 2
    Q[(2, 4)] = 2
    Q[(3, 5)] = 3
    Q[(4, 6)] = 3
    Q[(5, 6)] = 3
    Q[(6, 6)] = 3
    Q[(7, 5)] = 3
    Q[(8, 3)] = 2
    Q[(9, 2)] = 1
    Q[(9, 1)] = 1
    
    D = sorted(D, key=lambda p: p[0])
    Q = OrderedDict(sorted(Q.items(), key=lambda p: p[0][0]))
    P = GetClumpsPartition(D, Q)
    
    assert P[(1, 1)] == 1
    assert P[(1, 2)] == 1
    assert P[(1, 3)] == 1
    assert P[(1, 4)] == 1
    assert P[(2, 3)] == 2
    assert P[(2, 4)] == 2
    assert P[(3, 5)] == 3
    assert P[(4, 6)] == 3
    assert P[(5, 6)] == 3
    assert P[(6, 6)] == 3
    assert P[(7, 5)] == 3
    assert P[(8, 3)] == 4
    assert P[(9, 2)] == 5
    assert P[(9, 1)] == 5
    
    #Spinellis OpenMIC data
    D = [(0, 0), (1, 1), (3, 2), (2, 1), (5, 0), (4, 3), (6, 4)]
    D = sorted(D, key=lambda p: p[0])
    Q = EquipartitionYAxis(D, y=3)
    P = GetClumpsPartition(D, Q)
    clumps_partition_groups = GroupPartitionsPoints(P)
    assert (0, 0) in clumps_partition_groups[1]
    assert (1, 1) and (2, 1) in clumps_partition_groups[2]
    assert (3, 2) and (4, 3) in clumps_partition_groups[3]
    assert (5, 0) in clumps_partition_groups[4]
    assert (6, 4) in clumps_partition_groups[5]

def test_H():
    assert H(P=[0.25, 0.25, 0.25, 0.25]) == 2
    
    #OpenMIC test case
    assert H(P=[1./8, 1./4, 1./8, 1./2]) == 7./4

def test_getPartitionIndices():
    Q = EquipartitionYAxis(D, y=3)    
    P = GetClumpsPartition(D, OrderedDict(sorted(Q.items(), key=lambda p: p[0][0])))
    x_p = GetPartitionIndices(P, D, axis='x')
    y_p = GetPartitionIndices(Q, D, axis='y')    
    
    # Tested with visual inspection: See Albanese
    # http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf#page=3
    # visualize_partition(D, x_p, y_p)
    
    
test_EquipartitionYAxis()
test_GetClumpsPartition()
test_H()
test_getPartitionIndices()

Q = EquipartitionYAxis(D, y=3)
Q = OrderedDict(sorted(Q.items(), key=lambda p: p[0][0]))
P = GetClumpsPartition(D, Q)

