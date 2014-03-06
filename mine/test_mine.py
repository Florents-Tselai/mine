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
        x1 = np.array([1, 1, 1, 1, 2, 2, 3, 4, 5, 6 ,7, 8, 9, 9])
        y1 = np.array([1, 2, 3, 4, 3, 4, 5, 6, 6, 6, 5, 3, 2, 1])
        self.mine1 = MINE(x1,y1)

        
        # OpenMIC data
        x2 = np.array([0, 1, 3, 2, 5, 4, 6])
        y2 = np.array([0, 1, 2, 1, 0, 3, 4])
        self.mine2 = MINE(x2,y2)

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

    def test_get_clumps_partition(self):
    
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, 3)

        p1, ordinals1 = self.mine1.get_clumps_partition(q1)
        
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

        assert ordinals1 == [-1, 3, 5, 10, 11, 13]
        
        q2 = self.mine2.equipartition_y_axis(self.mine2.Dy, 3)
        p2, ordinals2 = self.mine2.get_clumps_partition(q2)
        assert p2[(0,0)] == 0

        assert p2[(1,1)] == 1
        assert p2[(2,1)] == 1

        assert p2[(3,2)] == 2
        assert p2[(4,3)] == 2

        assert p2[(5,0)] == 3

        assert p2[(6,4)] == 4

        assert ordinals2 == [-1, 0, 2, 4, 5, 6]

    def test_get_all_size_2_partition(self):
        ordinals = [2, 7, 9, 11, 14, 15]
        assert len(list(get_all_size_2_partition(ordinals))) == 14

    def test_HP(self):
        assert HP(np.array([2, 7, 9, 11, 14, 15])) == entropy([5./13, 2./13, 2./13, 3./13, 1./13])
        assert HP(np.array([-1, 0, 1])) == entropy([1./2, 1./2])

    def test_HQ(self):
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, y=3)
        assert HQ(q1) == entropy([4./14, 5./14, 5./14])

        q2 = self.mine2.equipartition_y_axis(self.mine2.Dy, y=3)
        assert HQ(q2) == entropy([2./7, 2./7, 3./7])

    def test_get_map_from_ordinals(self):
        assert self.mine1.get_map_from_ordinals(np.array([-1, 3, 8, 13]), axis='y') == self.mine1.equipartition_y_axis(self.mine1.Dy, 3)

        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, 3)

        p1, ordinals1 = self.mine1.get_clumps_partition(q1)
        assert p1 == self.mine1.get_map_from_ordinals(ordinals1, 'x')

    def test_get_grid_histogram(self):
        q1 = self.mine1.equipartition_y_axis(self.mine1.Dy, 3)
        p1, ordinals1 = self.mine1.get_clumps_partition(q1)
        #plot_partitions(p1, q1)
        #plt.show()
        #Inspect visually as well
        assert_array_equal(self.mine1.get_grid_histogram(ordinals1, q1), np.array([0,0,5,0,0,2,2,0,1,0,2,0,0,0,2]))

    def test_I(self):
        pass

    def test_MI(self):
        pass

    def test_number_of_points_in_partition(self):
        assert number_of_points_in_partition([2, 7, 9, 11, 14, 15]) == 13
        assert number_of_points_in_partition([-1, 2, 7, 9, 11, 14, 15]) == 16
        assert number_of_points_in_partition([-1, 7, 15]) == 16
        assert number_of_points_in_partition([-1, 0, 2, 5, 6]) == 7
        assert number_of_points_in_partition([-1,0,1]) == 2
        assert number_of_points_in_partition(([-1,1,2])) == 3
        assert number_of_points_in_partition([-1, 5, 6]) == 7
        assert number_of_points_in_partition([-1,0]) == 1
        assert number_of_points_in_partition([-1,0,1]) == 2

    def test_get_partition_histogram(self):
        assert_array_equal(get_partition_histogram(np.array([2, 7, 9, 11, 14, 15])), np.array([5 ,2, 2, 3, 1]))
        assert_array_equal(get_partition_histogram(np.array([-1, 2, 7, 9, 11, 14, 15])), np.array([3, 5, 2, 2, 3, 1]))
        assert_array_equal(get_partition_histogram(np.array([-1,7,15])), np.array([8, 8]))
        assert_array_equal(get_partition_histogram(np.array([-1, 0, 2, 5, 6])), np.array([1, 2, 3, 1]))
        assert_array_equal(get_partition_histogram(np.array([-1,0,1])), np.array([1, 1]))
        assert_array_equal(get_partition_histogram(np.array([-1,1,2])), np.array([2, 1]))
        assert_array_equal(get_partition_histogram(np.array([-1, 5, 6])), np.array([6, 1]))

    def test_optimize_x_axis(self):
        x = np.arange(30)
        y = x**2
        m = MINE(x,y)
        q = m.equipartition_y_axis(m.Dy, 5)
        m.optimize_x_axis(m.Dx, q, 4)

    def test_approx_char_matrix(self):
        x = np.linspace(0, 1, 100)
        y = np.sin(10 * np.pi * x) + x
        y +=np.random.uniform(-1, 1, x.shape[0])
        m = MINE(x, y)
        b = pow(len(x), 0.6)
        M = m.approx_char_matrix(m.D, b)
        print np.max(M[~np.isnan(M)])

    def test_comparison(self):
        from minepy import MINE


if __name__ == '__main__':
    unittest.main()
