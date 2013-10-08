# Copyright 2013-2014 Florents Tselai
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import OrderedDict

from mine import *


def test_EquipartitionYAxis():
    # Albanese et. al. data
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
    
    Q = EquipartitionYAxis(D, y=3)
    
    assert Q[(1, 1)] == 0
    assert Q[(1, 2)] == 0
    assert Q[(9, 2)] == 0
    assert Q[(9, 1)] == 0
    
    assert Q[(1, 3)] == 1
    assert Q[(1, 4)] == 1
    assert Q[(2, 3)] == 1
    assert Q[(2, 4)] == 1
    assert Q[(8, 3)] == 1
    
    assert Q[(3, 5)] == 2
    assert Q[(4, 6)] == 2
    assert Q[(5, 6)] == 2
    assert Q[(6, 6)] == 2
    assert Q[(7, 5)] == 2
    
    # Spinellis OpenMIC data
    D = [(0, 0), (1, 1), (3, 2), (2, 1), (5, 0), (4, 3), (6, 4)]
    
    '''
    -------------
    3         x
    2       x
    -------------
    1   x x
    -------------
    0 x         x
    -------------
      0 1 2 3 4 5
     '''
    Q = EquipartitionYAxis(D, y=3)
    
    assert Q[(0, 0)] == 0
    assert Q[(5, 0)] == 0
    
    assert Q[(1, 1)] == 1
    assert Q[(2, 1)] == 1
    
    assert Q[(3, 2)] == 2
    assert Q[(4, 3)] == 2
    
def test_GetClumpsPartition():
    
    # Albanese et.al. data
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]

    Q = EquipartitionYAxis(D, 3)
    
    P = GetClumpsPartition(D, Q)
    
    assert P[(1, 1)] == 0
    assert P[(1, 2)] == 0
    assert P[(1, 3)] == 0
    assert P[(1, 4)] == 0
    
    assert P[(2, 3)] == 1
    assert P[(2, 4)] == 1
    
    assert P[(3, 5)] == 2
    assert P[(4, 6)] == 2
    assert P[(5, 6)] == 2
    assert P[(6, 6)] == 2
    assert P[(7, 5)] == 2
    
    assert P[(8, 3)] == 3
    
    assert P[(9, 2)] == 4
    assert P[(9, 1)] == 4
    
    # Spinellis OpenMIC data
    D = [(0, 0), (1, 1), (3, 2), (2, 1), (5, 0), (4, 3), (6, 4)]
    Q = EquipartitionYAxis(D, y=3)
    P = GetClumpsPartition(D, Q)
    clumps_partition_groups = GroupPointsByPartition(P)
    assert (0, 0) in clumps_partition_groups[0]
    assert (1, 1) and (2, 1) in clumps_partition_groups[1]
    assert (3, 2) and (4, 3) in clumps_partition_groups[2]
    assert (5, 0) in clumps_partition_groups[3]
    assert (6, 4) in clumps_partition_groups[4]

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
     
def test_visualize():
    # Albanese dataset
    # http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf#page=3
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
    Q = EquipartitionYAxis(D, y=3)
    P = GetClumpsPartition(D, Q)
    
    visualize(P, Q)

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
    
    grid_matrix = GetGridMatrix(P, Q)
    assert grid_matrix[0][0] == 0
    assert grid_matrix[0][1] == 0
    assert grid_matrix[0][2] == 2
    assert grid_matrix[1][1] == 2
    assert grid_matrix[2][0] == 1
    assert grid_matrix[2][2] == 1

def test_GetOrdinals():
    # Points sorted by increasing x-value
    D = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 0), (6, 4)]
    
    # Given the following x-axis partition
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
     
    P1 = {
          (0, 0): 0,
          (1, 1): 1,
          (2, 1): 1,
          (3, 2): 2,
          (4, 3): 2,
          (5, 0): 2,
          (6, 4): 3
         }
    
    partition_size = len(set(P1.values()))
    assert partition_size == 4
    
    expected_ordinals = [-1, 0, 2, 5, 6]
    assert len(GetPartitionOrdinals(D, P1)) == partition_size + 1
    assert GetPartitionOrdinals(D, P1) == expected_ordinals
    
    # Another test
    """
      4            |x
      3         x  |
      2       x    |
         *----+---+-----+-
      1   x x      |
         *----+---+-----+-
      0 x         x|
        0 1 2 3 4 5 6
                   |
     """
    P2 = {
          (0, 0): 0,
          (1, 1): 0,
          (2, 1): 0,
          (3, 2): 0,
          (4, 3): 0,
          (5, 0): 0,
          (6, 4): 1
         }
    partition_size = len(set(P2.values()))
    assert partition_size == 2
    
    expected_ordinals = [-1, 5, 6]
    assert len(GetPartitionOrdinals(D, P2)) == partition_size + 1
    assert GetPartitionOrdinals(D, P2) == expected_ordinals
   
def test_GetPartitionFromOrdinals():
    # Points sorted by increasing x-value
    D = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 0), (6, 4)]
    
    # Given the following x-axis points ordinals
    ordinals = [-1, 0, 2, 5, 6]
    # We should get the following x-axis partition
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
    # Which is described by this dictionary
    P = {(0, 0): 0,
         (1, 1): 1,
         (2, 1): 1,
         (3, 2): 2,
         (4, 3): 2,
         (5, 0): 2,
         (6, 4): 3
         }
    assert GetPartitionFromOrdinals(D, ordinals) == P
    
    # Now test with two ordinals
    
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
     
    ordinals = [-1, 2, 5, 6]
    P = {(0, 0): 0,
         (1, 1): 0,
         (2, 1): 0,
         (3, 2): 1,
         (4, 3): 1,
         (5, 0): 1,
         (6, 4): 2
         }
    assert GetPartitionFromOrdinals(D, ordinals) == P
    
def test_OptimizeXAxis():
    pass

def test_ApproxCharacteristicMatrix():
    x = range(100)
    y = [a for a in x]
    D = zip(x, y)
    n = len(D)
    B = pow(n, 0.6)
    c = 15
    M = ApproxCharacteristicMatrix(D, B, c)
    print M
    
    

# Run all tests

def test_GetPartitionHistogram():
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
    D = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 0), (6, 4)]

    ordinals = [-1, 2, 5, 6]
    assignments, probs = GetPartitionHistogram(D, ordinals)
    assert list(assignments) == [3,3,1]
    assert sum(probs) == 1.0
    
    ############ Another Test Case ##############
    
    
     # Points sorted by increasing x-value
    D = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 0), (6, 4)]
    
    # Given the following x-axis points ordinals
    ordinals = [-1, 0, 2, 5, 6]
    # We should get the following x-axis partition
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
    assignments, probs = GetPartitionHistogram(D, ordinals)
    assert list(assignments) == [1,2,3,1]
    assert sum(probs) == 1.0
    

#test_EquipartitionYAxis()
#test_GetClumpsPartition()
#test_H()
# test_visualize()
#test_GetGridMatrix()
#test_GetOrdinals()
#test_GetPartitionFromOrdinals()
#test_OptimizeXAxis()
#test_ApproxCharacteristicMatrix()
test_GetPartitionHistogram()