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

import numpy as np
from numpy.testing import assert_array_equal
import matplotlib.pyplot as plt

from mine import *


class mine__test(unittest.TestCase):
    def setUp(self):
        # Albanese et. al. data
        x1 = np.array([1, 1, 1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9])
        y1 = np.array([1, 2, 3, 4, 3, 4, 5, 6, 6, 6, 5, 3, 2, 1])
        self.mine1 = MINE(x1, y1)


        # OpenMIC data
        x2 = np.array([0, 1, 3, 2, 5, 4, 6])
        y2 = np.array([0, 1, 2, 1, 0, 3, 4])
        self.mine2 = MINE(x2, y2)

    def test_equipartition_y_axis(self):
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, 3)

        assert q1[(1, 1)] == 0
        assert q1[(1, 2)] == 0
        assert q1[(9, 2)] == 0
        assert q1[(9, 1)] == 0

        assert q1[(1, 3)] == 1
        assert q1[(1, 4)] == 1
        assert q1[(2, 3)] == 1
        assert q1[(2, 4)] == 1
        assert q1[(8, 3)] == 1

        assert q1[(3, 5)] == 2
        assert q1[(4, 6)] == 2
        assert q1[(5, 6)] == 2
        assert q1[(6, 6)] == 2
        assert q1[(7, 5)] == 2


        '''
        4            x
        --------------
        3         x
        2       x
        ---------------
        1   x x
        ---------------
        0 x         x
        ---------------
          0 1 2 3 4 5 6
         '''
        q2 = self.mine2.equipartition_y_axis(self.mine2.Dy, y=3)
        assert q2[(0, 0)] == 0
        assert q2[(5, 0)] == 0

        assert q2[(1, 1)] == 1
        assert q2[(2, 1)] == 1

        assert q2[(3, 2)] == 2
        assert q2[(4, 3)] == 2
        assert q2[(6, 4)] == 2

        #assert len(q2.ordinals) == 4

    def test_get_clumps_partition(self):
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, 3)

        p1 = self.mine1.get_clumps_partition(q1)


        assert p1[(1, 1)] == 0
        assert p1[(1, 2)] == 0
        assert p1[(1, 3)] == 0
        assert p1[(1, 4)] == 0

        assert p1[(2, 3)] == 1
        assert p1[(2, 4)] == 1

        assert p1[(3, 5)] == 2
        assert p1[(4, 6)] == 2
        assert p1[(5, 6)] == 2
        assert p1[(6, 6)] == 2
        assert p1[(7, 5)] == 2

        assert p1[(8, 3)] == 3

        assert p1[(9, 2)] == 4
        assert p1[(9, 1)] == 4


        q2 = self.mine2.equipartition_y_axis(self.mine2.Dy, 3)
        p2 = self.mine2.get_clumps_partition(q2)

        assert p2[(0, 0)] == 0

        assert p2[(1, 1)] == 1
        assert p2[(2, 1)] == 1

        assert p2[(3, 2)] == 2
        assert p2[(4, 3)] == 2

        assert p2[(5, 0)] == 3

        assert p2[(6, 4)] == 4

    def test_get_super_clumps_partition(self):
        x = np.arange(0,100)
        y = x+x**2
        m = MINE(x,y)
        q = m.equipartition_y_axis(m.Dy, 20)
        clumps = m.get_clumps_partition(q)
        c_clumps = m.get_c(clumps)
        k_clumps = len(c_clumps) - 1
        k_hat = 7
        super_clumps = m.get_super_clumps_partition(q, k_hat)
        c_super_clumps = m.get_c(super_clumps)
        k_super_clumps = len(c_super_clumps)-1
        assert k_super_clumps <= k_hat
        assert clumps.keys() == super_clumps.keys()

    def test_hp(self):
        assert self.mine1.hp([2, 7, 9, 11, 13]) == entropy([5./11, 2./11, 2./11, 2./11])
        assert self.mine1.hp([-1, 0, 1]) == entropy([1./2, 1./2])

    def test_hq(self):
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, y=3)
        assert self.mine1.hq(q1) == entropy([4. / 14, 5. / 14, 5. / 14])

        q2 = self.mine2.equipartition_y_axis(self.mine2.Dy, y=3)
        assert self.mine2.hq(q2) == entropy([2. / 7, 2. / 7, 3. / 7])

    def test_extend(self):
        x = np.arange(100)
        y = x**2 -2*x + 2
        m = MINE(x,y)

        assert_array_equal(m.extend(np.array([-1,5,7,8,10,15]), 6), np.array([-1,5,6,7,8,10,15]))

        assert_array_equal(m.extend(np.array([-1,5,7,8,10,15]), 7), np.array([-1,5,7,8,10,15]))

        assert_array_equal(m.extend(np.array([-1,5,7,8,10,15]), 17), np.array([-1,5,7,8,10,15,17]))

    def test_c(self):
        x = np.arange(10)
        y = x
        m = MINE(x,y)
        p = {(0, 0):0,
             (1, 1):0,
             (2, 2):0,
             (3, 3):1,
             (4, 4):1,
             (5, 5):2,
             (6, 6):3,
             (7, 7):3,
             (8, 8):3,
             (9, 9):4
            }
        assert m.get_c(p) == [-1, 2, 4, 5, 8, 9]

        p = {
             (3, 3):0,
             (4, 4):0,
             (5, 5):1,
             (6, 6):2,
             (7, 7):2,
             (8, 8):2,
             (9, 9):3
            }

        assert m.get_c(p) == [2, 4, 5, 8, 9]

    def test_get_grid_histogram(self):
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, 3)
        p1 = self.mine1.get_clumps_partition(q1)

        #plot_partitions(p1, q1)
        #plt.show()
        #Inspect visually as well
        assert_array_equal(self.mine1.get_grid_histogram(p1, q1)
                           ,np.array([0, 0, 5, 0, 0, 2, 2, 0, 1, 0, 2, 0, 0, 0, 2]))

    def test_optimize_x_axis(self):

        x = np.arange(1000)
        y = x**2 -2*x + 2
        m = MINE(x,y)
        q = m.equipartition_y_axis(m.Dy, 5)
        x_size = 10
        opt = m.optimize_x_axis(m.Dx,q, x_size)
        assert len(opt) == x_size - 1




if __name__ == '__main__':
    unittest.main()
