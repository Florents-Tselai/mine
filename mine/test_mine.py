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
        #Albanese et. al. data
        self.D1 = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 6), (7, 5), (8, 3), (9, 2), (9, 1)]
        self.D1y = sort_D_increasing_by(self.D1, 'y')
        self.D1x = sort_D_increasing_by(self.D1, 'x')
        
        #OpenMIC data
        self.D2 = [(0, 0), (1, 1), (3, 2), (2, 1), (5, 0), (4, 3), (6, 4)]
        self.D2y = sort_D_increasing_by(self.D2, 'y')
        self.D2x = sort_D_increasing_by(self.D2, 'x')
          
    def test_EquipartitionYAxis(self):
        
        Q = EquipartitionYAxis(self.D1y, y=3)
        
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
        Q = EquipartitionYAxis(self.D2y, y=3)
        
        assert Q[(0, 0)] == 0
        assert Q[(5, 0)] == 0
        
        assert Q[(1, 1)] == 1
        assert Q[(2, 1)] == 1
        
        assert Q[(3, 2)] == 2
        assert Q[(4, 3)] == 2
        
    def test_GetClumpsPartition(self):
        
    
        Q = EquipartitionYAxis(self.D1y, 3)
        
        P = GetClumpsPartition(self.D1x, Q)
        
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
        
        Q = EquipartitionYAxis(self.D2y, y=3)
        P = GetClumpsPartition(self.D2x, Q)
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
          """
          4  |   |     |x|
          3  |   |  x  | |
          2  |   |x    | |
          1  |x x|     | |
          0 x|   |    x| |
            0|1 2|3 4 5|6|
          """
          assert HP(np.array([-1, 0, 2, 5, 6])) == H(np.array([1. / 7, 2. / 7 , 3. / 7, 1. / 7]))
          assert HP(np.array([0, 2, 5])) == H(np.array([2. / 5, 3. / 5]))
          assert HP(np.array([-1, 0])) == H(np.array([1. / 1]))
          assert HP(np.array([-1, 5, 6])) == H(np.array([6. / 7 , 1. / 7]))
        
    def test_HQ(self):
        """
          4             x
          3         x   
          2       x     
            *----+---+-----+-
          1   x         
             *----+---+-----+-
          0 x        x
            0 1 2 3 4 5 6
                    
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
                ) == H(np.array([2. / 6, 1. / 6, 3. / 6]))
                
        """
          4             x
          3         x   
          2       x     
             *----+---+-----+-
          1   x x       
          0 x         x 
            0 1 2 3 4 5 6
                      
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
                ) == H(np.array([4. / 7, 3. / 7]))
        
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
        
        self.assertEqual(HPQ(P_ordinals, Q_map), H(np.array([0. / 6, 0. / 6, 2. / 6, 0. / 6, 2. / 6, 0. / 6, 1. / 6, 0. / 6, 1./6])))
    

if __name__ == '__main__':
    unittest.main()