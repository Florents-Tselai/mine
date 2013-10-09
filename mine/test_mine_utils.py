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
    
def test_visualize():
    # Albanese dataset
    # http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf#page=3
    D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
    Q = EquipartitionYAxis(D, y=3)
    P = GetClumpsPartition(D, Q)
    
    visualize(P, Q)
    
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
    