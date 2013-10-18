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

import unittest
from collections import OrderedDict

from mine import *

class mine__test(unittest.TestCase):

    def setUp(self):
        pass
    
    
    def test_EquipartitionYAxis(self):
        # Albanese et. al. data
        D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
        
        Dy  = sort_D_increasing_by(D, 'y')
        Q = EquipartitionYAxis(sort_D_increasing_by(D, 'y'), y=3)
        
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
        Q = EquipartitionYAxis(sort_D_increasing_by(D,increasing_by='y'), y=3)
        
        assert Q[(0, 0)] == 0
        assert Q[(5, 0)] == 0
        
        assert Q[(1, 1)] == 1
        assert Q[(2, 1)] == 1
        
        assert Q[(3, 2)] == 2
        assert Q[(4, 3)] == 2
        
    def test_GetClumpsPartition(self):
        
        # Albanese et.al. data
        D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
    
        Q = EquipartitionYAxis(sort_D_increasing_by(D,increasing_by='y'), 3)
        
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
        Q = EquipartitionYAxis(sort_D_increasing_by(D,increasing_by='y'), y=3)
        P = GetClumpsPartition(sort_D_increasing_by(D,increasing_by='x'), Q)
        clumps_partition_groups = GroupPointsByPartition(P)
        assert (0, 0) in clumps_partition_groups[0]
        assert (1, 1) and (2, 1) in clumps_partition_groups[1]
        assert (3, 2) and (4, 3) in clumps_partition_groups[2]
        assert (5, 0) in clumps_partition_groups[3]
        assert (6, 4) in clumps_partition_groups[4]
    
    def test_OptimizeXAxis(self):
        pass
    
    def test_ApproxCharacteristicMatrix(self):
        x = range(100)
        y = [a ^ 2 + a ^ 3 + 2 * a for a in x]
        D = zip(x, y)
        n = len(D)
        B = pow(n, 0.6)
        c = 15
        M = ApproxCharacteristicMatrix(D, B, c)
        print M
        
    def test_HP(self):
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
        assert HP(Dx, P1_ordinals) == entropy([1. / 7, 2. / 7 , 3. / 7, 1. / 7])
        
        assert HP(Dx, [0, 2, 5]) == entropy([2. / 5, 3. / 5])
        
        assert HP(Dx, [-1, 0]) == entropy([1. / 1])
        
        assert HP(Dx, [-1, 5, 6]) == entropy([6. / 7 , 1. / 7])
        
    def test_HQ(self):
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
                 (0, 0): 0,
                 (5, 0): 0,
                 (1, 1): 1,
                 (3, 2): 2,
                 (4, 3): 2,
                 (6, 4): 2
                }
                ) == entropy([2. / 6, 1. / 6, 3. / 6])
                
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
                 (0, 0): 0,
                 (5, 0): 0,
                 (1, 1): 0,
                 (2, 1): 0,
                 (3, 2): 1,
                 (4, 3): 1,
                 (6, 4): 1
                }
                ) == entropy([4. / 7, 3. / 7])
        
    def test_HPQ(self):
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
        P_ordinals = [-1, 0, 2, 5]
        Q_map = {
                 (0, 0): 0,
                 (5, 0): 0,
                 
                 (1, 1): 1,
                 (2, 1): 1,
                 
                 (3, 2): 2,
                 (4, 3): 2,
                 (6, 4): 2
                 }
        
        assert HPQ(P_ordinals, Q_map) == entropy([1. / 6, 0. / 6, 0. / 6, 0. / 6, 2. / 6, 1. / 6, 0. / 6, 2. / 6])
     
    def test_GetPartitionOrdinalsFromMap(self):
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
        
        D = [(0, 0), (10, 10), (20, 20), (30, 30), (40, 40), (50, 50), (60, 60), (70, 70), (80, 80), (90, 90)]
        
        Q = {
             (0, 0): 0,
             (10, 10): 0,
             (20, 20): 0,
             (30, 30): 1,
             (40, 40): 1,
             (50, 50): 1,
             (60, 60): 2,
             (70, 70): 2,
             (80, 80): 3,
             (90, 90):4
             }
        
        expected_ordinals = [-1, 2, 5, 7, 8, 9]
        assert GetPartitionOrdinalsFromMap(D, Q) == expected_ordinals
        
    def test_GetParitionMapFromOrdinals(self):
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
        
    def test_GetPartitionHistogram(self):
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
        assert list(assignments) == [3, 3, 1]
        
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
        assert list(assignments) == [1, 2, 3, 1]
        
        # Other test case
        D = [(0, 0), (10, 10), (20, 20), (30, 30), (40, 40), (50, 50), (60, 60), (70, 70), (80, 80), (90, 90)]
        
        Q = {
             (0, 0): 0,
             (10, 10): 0,
             (20, 20): 0,
             (30, 30): 1,
             (40, 40): 1,
             (50, 50): 1,
             (60, 60): 2,
             (70, 70): 2,
             (80, 80): 3,
             (90, 90):4
             }
        
        assignments = GetPartitionHistogram(D, GetPartitionOrdinalsFromMap(D, Q, 'y'), 'y')
        assert assignments == [3, 3, 2, 1, 1]
     
    def test_visualize(self):
        # Albanese dataset
        # http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf#page=3
        D = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
        Q = EquipartitionYAxis(sort_D_increasing_by(D,increasing_by='y'), y=3)
        P = GetClumpsPartition(D, Q)
        
        #visualize(P, Q)   
     
    def test_GetGridHistogram(self):
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
        Dy = sort_D_increasing_by(D, increasing_by='y')
        Dx = sort_D_increasing_by(D, increasing_by='x')
        
        Q1 = {
             (0, 0): 0,
             (5, 0): 0,
             (1, 1): 1,
             (2, 1): 1,
             (3, 2): 2,
             (4, 3): 2,
             (6, 4):2
             }
        
        P1 = [-1, 0, 2, 5, 6]
        
        assert GetGridHistogram(Q1, P1) == [1, 0, 0, 0, 2, 0, 1, 0, 2, 0, 0, 1]
        
        P2 = [-1, 0, 5, 6]
        assert GetGridHistogram(Q1, P2) == [1, 0, 0, 1, 2, 2, 0, 0, 1]
    
    def test_m(self):
        assert m([-1, 2, 3, 6]) == 7
        assert m([2,5,9]) == 7
        assert m([-1,6,7]) == 8
        


if __name__ == '__main__':
    unittest.main()