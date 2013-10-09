from collections import defaultdict, Mapping
from mine_utils import *
import numpy as np

def test_GetPartitionOrdinalsFromMap():
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
    assert len(GetPartitionOrdinalsFromMap(D, P1)) == partition_size + 1
    assert GetPartitionOrdinalsFromMap(D, P1) == expected_ordinals
    
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
    assert len(GetPartitionOrdinalsFromMap(D, P2)) == partition_size + 1
    assert GetPartitionOrdinalsFromMap(D, P2) == expected_ordinals
    
    D = [(0,0), (10,10), (20,20), (30,30), (40,40), (50,50), (60,60), (70,70), (80,80), (90,90)]
    
    Q = {
         (0,0): 0, 
         (10,10): 0, 
         (20,20): 0, 
         (30,30): 1, 
         (40,40): 1, 
         (50,50): 1, 
         (60,60): 2, 
         (70,70): 2, 
         (80,80): 3, 
         (90,90):4
         }
    
    expected_ordinals = [-1, 2, 5, 7, 8, 9]
    assert GetPartitionOrdinalsFromMap(D, Q) == expected_ordinals
    
def test_GetParitionMapFromOrdinals():
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
    assert GetPartitionMapFromOrdinals(D, ordinals) == P
    
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
    assert GetPartitionMapFromOrdinals(D, ordinals) == P
    
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
    assignments = GetPartitionHistogram(D, ordinals)
    assert list(assignments) == [3,3,1]
    
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
    assignments = GetPartitionHistogram(D, ordinals)
    assert list(assignments) == [1,2,3,1]
    
    #Other test case
    D = [(0,0), (10,10), (20,20), (30,30), (40,40), (50,50), (60,60), (70,70), (80,80), (90,90)]
    
    Q = {
         (0,0): 0, 
         (10,10): 0, 
         (20,20): 0, 
         (30,30): 1, 
         (40,40): 1, 
         (50,50): 1, 
         (60,60): 2, 
         (70,70): 2, 
         (80,80): 3, 
         (90,90):4
         }
    
    assignments = GetPartitionHistogram(D, GetPartitionOrdinalsFromMap(D, Q, 'y'), 'y')
    assert assignments == [3,3,2,1,1]
 
def test_visualize():
    # Albanese dataset
    # http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf#page=3
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
    Q = EquipartitionYAxis(D, y=3)
    P = GetClumpsPartition(D, Q)
    
    visualize(P, Q)   
 
def test_GetGridHistogram():
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
    D = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 0), (6, 4)]
    Dy = sort_D_increasing_by(D, 'y')
    Dx = sort_D_increasing_by(D, 'x')
    
    Q1 = {
         (0,0): 0,
         (5,0): 0,
         (1,1): 1,
         (2,1): 1,
         (3,2): 2,
         (4,3): 2,
         (6,4):2
         }
    
    P1 = {
         (3, 2): 1, 
         (5, 0): 1, 
         (2, 1): 0, 
         (4, 3): 1, 
         (1, 1): 0
         }
    
    assert GetGridHistogram(Q1,P1) == [0,2,2,0,0,1]
    
    P2 = {
         (3, 2): 2, 
         (5, 0): 2, 
         (2, 1): 1,
         (1, 1): 1, 
         (4, 3): 2, 
         (0, 0): 0,
         (6, 4): 3
         }
    
    assert GetGridHistogram(Q1,P2) == [0,0,2,1,0,2,0,0,1,0,1,0]
    
     
test_GetPartitionOrdinalsFromMap()
test_GetParitionMapFromOrdinals()
test_GetPartitionHistogram()
test_GetGridHistogram()