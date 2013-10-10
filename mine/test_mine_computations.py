from mine_computations import *
from mine_utils import *

def test_HP():
    D = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 0), (6, 4)]
    Dx = sort_D_increasing_by(D, 'x')
    
    P1_ordinals = [-1, 0, 2, 5, 6]
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
    assert HP(Dx, P1_ordinals) == entropy([1./7, 2./7 , 3./7,1./7])
    
    assert HP(Dx, [0,2,5]) == entropy([2./5, 3./5])
    
    assert HP(Dx, [-1,0]) == entropy([1./1])
    
    assert HP(Dx, [-1, 5, 6]) == entropy([6./7 , 1./7])
    
def test_HQ():
    """
      4  |   |     |x
      3  |   |  x  |
      2  |   |x    |
         *----+---+-----+-
      1  |x  |     |
         *----+---+-----+-
      0 x|   |    x|
        0|1 2|3 4 5 6
         |   |     |
     """
     
    assert HQ(
            {
             (0,0): 0,
             (5,0): 0,
             (1,1): 1,
             (3,2): 2,
             (4,3): 2,
             (6,4): 2
            }
            ) == entropy([2./6, 1./6, 3./6])
            
    """
      4  |   |     |x
      3  |   |  x  |
      2  |   |x    |
         *----+---+-----+-
      1  |x x|     |
      0 x|   |    x|
        0|1 2|3 4 5 6
         |   |     |
     """
    assert HQ(
            {
             (0,0): 0,
             (5,0): 0,
             (1,1): 0,
             (2,1): 0,
             (3,2): 1,
             (4,3): 1,
             (6,4): 1
            }
            ) == entropy([4./7, 3./7])
    
def test_HPQ():
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
    P_ordinals = [0,2,5]
    Q_map = {
             (0,0): 0,
             (5,0): 0,
             (1,1): 1,
             (2,1): 1,
             (3,2): 2,
             (4,3): 2,
             (6,4): 2
             }
    Dx = sort_D_increasing_by(Q_map.keys(), 'x')
    GetPartitionMapFromOrdinals(Dx, P_ordinals, axis='x')
    assert HPQ(P_ordinals, Q_map) == entropy([0./5, 2./5, 2./5, 0./5, 1./5])

test_HP()
test_HQ()
#test_HPQ()