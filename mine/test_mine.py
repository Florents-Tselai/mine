from collections import OrderedDict
from mine import *

def test_EquipartitionYAxis():
    # Albanese data
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
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
    
    # Spinellis OpenMIC data
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
    
    # Spinellis OpenMIC data
    D = [(0, 0), (1, 1), (3, 2), (2, 1), (5, 0), (4, 3), (6, 4)]
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
    P = {(0, 0): 1,
         (1, 1): 2,
         (2, 1): 2,
         (3, 2): 3,
         (4, 3): 3,
         (5, 0): 3,
         (6, 4): 4
         }
    
    Q = {(0, 0): 1,
         (1, 1): 2,
         (2, 1): 2,
         (3, 2): 3,
         (4, 3): 3,
         (5, 0): 1,
         (6, 4): 3
         }
    
    # visualize_grid(P, Q)
    
    # Joint entropy computation
    assert (H(P, Q) == H([0 , 0   , 2. / 7, 1. / 7,
                         0 , 2. / 7 , 0. / 7, 0   ,
                        1. / 7, 0   , 1. / 7, 0   ]))
     
def test_visualize_grid():
    # Albanese dataset
    # http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf#page=3
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
    Q = EquipartitionYAxis(D, y=3)
    P = GetClumpsPartition(D, Q)
    
    # visualize_grid(P, Q)

def test_GetGridMatrix():
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
    P = {(0, 0): 1,
         (1, 1): 2,
         (2, 1): 2,
         (3, 2): 3,
         (4, 3): 3,
         (5, 0): 3,
         (6, 4): 4
         }
    
    Q = {(0, 0): 1,
         (1, 1): 2,
         (2, 1): 2,
         (3, 2): 3,
         (4, 3): 3,
         (5, 0): 1,
         (6, 4): 3
         }
    
    grid_matrix = GetGridMatrix(P, Q)
    assert grid_matrix[0][0] == 0
    assert grid_matrix[0][1] == 0
    assert grid_matrix[0][2] == 2
    assert grid_matrix[1][1] == 2
    assert grid_matrix[2][0] == 1
    assert grid_matrix[2][2] == 1

# Run tests

def test_GetOrdinals():
    #Points sorted by increasing x-value
    D = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 0), (6, 4)]
    
    #Given the following x-axis partition
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
     
    P = {(0, 0): 1,
         (1, 1): 2,
         (2, 1): 2,
         (3, 2): 3,
         (4, 3): 3,
         (5, 0): 3,
         (6, 4): 4
         }
    
    #The endpoints are (0,0), (2,1), (5,0) and the corresponding ordinals are:
    expected_ordinals = [0, 2, 5]
    assert GetPartitionOrdinals(D, P) ==  expected_ordinals
   
def test_GetPartitionFromOrdinals():
    #Points sorted by increasing x-value
    D = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 0), (6, 4)]
    
    #Given the following x-axis points ordinals
    ordinals = [0, 2, 5]
    #We should get the following x-axis partition
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
    #Which is described by this dictionary
    P = {(0, 0): 1,
         (1, 1): 2,
         (2, 1): 2,
         (3, 2): 3,
         (4, 3): 3,
         (5, 0): 3,
         (6, 4): 4
         }

    assert GetPartitionFromOrdinals(D, ordinals) == P
    
    #Now test with two ordinals
    
    """
      4      |     |x
      3      |  x  |
      2      |x    |
         *----+---+-----+-
      1  x x |     |
         *----+---+-----+-
      0 x    |    x|
        0 1 2|3 4 5 6
             |     |
     """
     
    ordinals = [2, 5]
    P = {(0, 0): 1,
         (1, 1): 1,
         (2, 1): 1,
         (3, 2): 2,
         (4, 3): 2,
         (5, 0): 2,
         (6, 4): 3
         }
    assert GetPartitionFromOrdinals(D, ordinals) == P
    
def test_OptimizeXAxis():
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
    Q = EquipartitionYAxis(D, y=3)
    I = OptimizeXAxis(D, Q, x=2, k_hat=1)
    print I
    

#Run all tests
test_EquipartitionYAxis()
test_GetClumpsPartition()
test_H()
test_visualize_grid()
test_GetGridMatrix()
test_GetOrdinals()
test_GetPartitionFromOrdinals()
test_OptimizeXAxis()