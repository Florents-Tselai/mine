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


    def test_hp(self):
        assert self.mine1.create_partition([2, 7, 9, 11, 13]).h() == entropy([5./11, 2./11, 2./11, 2./11])
        assert self.mine1.create_partition([-1, 0, 1]).h() == entropy([1./2, 1./2])

    def test_hq(self):
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, y=3)
        assert hq(q1) == entropy([4. / 14, 5. / 14, 5. / 14])

        q2 = self.mine2.equipartition_y_axis(self.mine2.Dy, y=3)
        assert hq(q2) == entropy([2. / 7, 2. / 7, 3. / 7])


    def test_get_grid_histogram(self):
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, 3)
        p1 = self.mine1.get_clumps_partition(q1)

        #plot_partitions(p1, q1)
        #plt.show()
        #Inspect visually as well
        assert_array_equal(p1.grid_histogram(q1)
                           ,np.array([0, 0, 5, 0, 0, 2, 2, 0, 1, 0, 2, 0, 0, 0, 2]))

    def test_number_of_points(self):
        assert self.mine1.create_partition([-1, 7, 13]).number_of_points()  == 8 + 6
        assert self.mine1.create_partition([-1, 0, 2, 5, 6]).number_of_points()  == 1 + 2 + 3 + 1
        assert self.mine1.create_partition([-1, 0, 1]).number_of_points()  == 1 + 1
        assert self.mine1.create_partition(([-1, 1, 2])).number_of_points() == 2 + 1
        assert self.mine1.create_partition([-1, 5, 6]).number_of_points()  == 6 + 1
        assert self.mine1.create_partition([-1, 0]).number_of_points()  == 1
        assert self.mine1.create_partition([-1, 0, 1]).number_of_points()  == 1 + 1
        assert self.mine1.create_partition([-1, 3, 5, 10, 11, 13]).number_of_points() == 4 + 2 + 5 + 1 + 2
        assert self.mine1.create_partition([2, 7, 9, 11, 13]).number_of_points() == 5 + 2 + 2 + 2


    def test_histogram(self):
        assert_array_equal(self.mine1.create_partition(np.array([2, 7, 9, 11, 13])).histogram(), np.array([5, 2, 2, 2]))
        assert_array_equal(self.mine1.create_partition(np.array([-1, 2, 7, 9, 11, 13])).histogram(), np.array([3, 5, 2, 2, 2]))
        assert_array_equal(self.mine1.create_partition(np.array([-1, 7, 13])).histogram(), np.array([8, 6]))
        assert_array_equal(self.mine1.create_partition(np.array([-1, 0, 2, 5, 6])).histogram(), np.array([1, 2, 3, 1]))
        assert_array_equal(self.mine1.create_partition(np.array([-1, 0, 1])).histogram(), np.array([1, 1]))
        assert_array_equal(self.mine1.create_partition(np.array([-1, 1, 2])).histogram(), np.array([2, 1]))
        assert_array_equal(self.mine1.create_partition(np.array([-1, 5, 6])).histogram(), (np.array([6, 1])))

    def test_get_super_clumps_partition(self):
        x = np.arange(100)
        y = x**2-3*x**3 + np.sqrt(x)
        m = MINE(x,y)
        q = m.equipartition_y_axis(m.Dy, 10)
        #clumps = m.get_clumps_partition(q)
        #superclumps = m.get_super_clumps_partition(q, 8)

        #assert len(superclumps) == 8
        #assert clumps.points() == superclumps.points()

if __name__ == '__main__':
    unittest.main()
