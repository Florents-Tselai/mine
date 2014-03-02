# Copyright 2014 Florents Tselai
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
import unittest

from mine import *


class mine__test(unittest.TestCase):

    def setUp(self):
        # Albanese et. al. data
        x1 = np.array([1, 1, 1, 1, 2, 2, 3, 4, 5, 6 ,7, 8, 9, 9])
        y1 = np.array([1, 2, 3, 4, 3, 4, 5, 6, 6, 6, 5, 3, 2, 1])
        self.mine1 = MINE(x1,y1)

        
        # OpenMIC data
        x2 = np.array([0, 1, 3, 2, 5, 4, 6])
        y2 = np.array([0, 1, 2, 1, 0, 3, 4])
        self.mine2 = MINE(x2,y2)

    def test_equipartition_y_axis(self):

        Q = self.mine1.get_points_assignments(self.mine1.equipartition_y_axis(y=3))

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
        Q = {(self.mine2.Dy[k][0], self.mine2.Dy[k][1]): v for k,v in self.mine2.equipartition_y_axis(y=3).iteritems()}
        
        assert Q[(0, 0)] == 0
        assert Q[(5, 0)] == 0
        
        assert Q[(1, 1)] == 1
        assert Q[(2, 1)] == 1
        
        assert Q[(3, 2)] == 2
        assert Q[(4, 3)] == 2
        
    def test_get_clumps_partition(self):
    
        q = self.mine1.equipartition_y_axis(3)

        p = self.mine1.get_clumps_partition(q)
        
        assert p[(1, 1)] == 0
        assert p[(1, 2)] == 0
        assert p[(1, 3)] == 0
        assert p[(1, 4)] == 0
        
        assert p[(2, 3)] == 1
        assert p[(2, 4)] == 1
        
        assert p[(3, 5)] == 2
        assert p[(4, 6)] == 2
        assert p[(5, 6)] == 2
        assert p[(6, 6)] == 2
        assert p[(7, 5)] == 2
        
        assert p[(8, 3)] == 3
        
        assert p[(9, 2)] == 4
        assert p[(9, 1)] == 4
        
        q = self.mine2.equipartition_y_axis(3)
        p = self.mine2.get_clumps_partition(q)
        clumps_partition_groups = GroupPointsByPartition(p)
        assert (0, 0) in clumps_partition_groups[0]
        assert (1, 1) and (2, 1) in clumps_partition_groups[1]
        assert (3, 2) and (4, 3) in clumps_partition_groups[2]
        assert (5, 0) in clumps_partition_groups[3]
        assert (6, 4) in clumps_partition_groups[4]
                          
    def test_mine(self):
        return
        x = np.array(range(1000))
        y = 4 * (x - 1. / 2) ** 2
        D = zip(x, y)
        n = len(D)
        B = pow(n, 0.6)
        c = 15
        M = ApproxCharacteristicMatrix(D, B, c=1)
        print mine(M, B, c)

if __name__ == '__main__':
    unittest.main()
