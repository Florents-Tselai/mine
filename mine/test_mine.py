from collections import OrderedDict
from mine import visualize_partition, EquipartitionYAxis, GetClumpsPartition

"""
Example from Albanese et. al. pg. 4
http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf

"""

D = [(1,1), (1,2), (1,3), (1,4), (2,3), (2,4), (3,5), (4,6), (5,6), (6,6), (7,5), (8,3), (9,2), (9,1)]
x_partition = [1.8, 2.2, 7.8, 8.2]
y_partition = [2.5, 4.6]

#visualize_partition(D, x_partition, y_partition)

def test_EquipartitionYAxis():
    Q = EquipartitionYAxis(D, y=3)
    
    assert Q[(1,1)] == 1
    assert Q[(1,2)] == 1
    assert Q[(1,3)] == 2
    assert Q[(1,4)] == 2
    assert Q[(2,3)] == 2
    assert Q[(2,4)] == 2
    assert Q[(3,5)] == 3
    assert Q[(4,6)] == 3
    assert Q[(5,6)] == 3
    assert Q[(6,6)] == 3
    assert Q[(7,5)] == 3
    assert Q[(8,3)] == 2
    assert Q[(9,2)] == 1
    assert Q[(9,1)] == 1    
    
def test_GetClumpsPartition():
    Q = {}
    Q[(1,1)] = 1
    Q[(1,2)] = 1
    Q[(1,3)] = 2
    Q[(1,4)] = 2
    Q[(2,3)] = 2
    Q[(2,4)] = 2
    Q[(3,5)] = 3
    Q[(4,6)] = 3
    Q[(5,6)] = 3
    Q[(6,6)] = 3
    Q[(7,5)] = 3
    Q[(8,3)] = 2
    Q[(9,2)] = 1
    Q[(9,1)] = 1
    
    Q = OrderedDict(sorted(Q.items(), key=lambda p: p[0][0]))
    P = GetClumpsPartition(D, Q)
    
    assert P[(1,1)] == 1
    assert P[(1,2)] == 1
    assert P[(1,3)] == 1
    assert P[(1,4)] == 1
    assert P[(2,3)] == 2
    assert P[(2,4)] == 2
    assert P[(3,5)] == 3
    assert P[(4,6)] == 3
    assert P[(5,6)] == 3
    assert P[(6,6)] == 3
    assert P[(7,5)] == 3
    assert P[(8,3)] == 4
    assert P[(9,2)] == 5
    assert P[(9,1)] == 5 

test_EquipartitionYAxis()
test_GetClumpsPartition()

Q = EquipartitionYAxis(D, y=3)
P = GetClumpsPartition(D,  OrderedDict(sorted(Q.items(), key=lambda p: p[0][0])))
for p in Q.iterkeys():
    print p, Q[p], P[p]
   