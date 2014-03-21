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

from __future__ import division
import bisect
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict, Counter
from copy import copy
from itertools import *
from math import floor
from matplotlib import cm
from numpy import vstack, lexsort, shape, where, log2, fliplr, float64
import matplotlib.pyplot as plt
import numpy as np

p_x, p_y = lambda p: p[0], lambda p: p[1]


class MINE:
    def __init__(self, x, y):
        self.D = np.core.records.fromarrays([x, y], names='x,y')
        self.D_orth = list(np.core.records.fromarrays([y, x], names='x,y'))
        self.Dx = [tuple(p) for p in self.D[self.D.argsort(order='x')]]
        self.Dy = [tuple(p) for p in self.D[self.D.argsort(order='y')]]
        self.n = len(self.D)

    def approx_max_mi(self, d, x, y):
        q = self.equipartition_y_axis(self.Dy, y)
        #return self.optimize_x_axis(self.Dx, q, x, k_hat)
        return self.optimize_x_axis(self.Dx, q, x)

    def approx_char_matrix(self, d, b):
        from math import ceil, floor

        s = int(floor(b / 2)) + 1

        I = np.zeros(shape=(s, s), dtype=float64)
        I_orth = np.zeros(shape=(s, s), dtype=float64)
        M = np.zeros(shape=(s, s), dtype=float64)

        for y in range(2, int(floor(b / 2)) + 1):
            x = int(floor(b / y))
            assert len(self.approx_max_mi(self.D, x, y)) == x - 1
            print x, y
            #print "==="
            #print len(I[2:x+1][y]), len(self.approx_max_mi(self.D, x, y))
            I[2:x + 1][y] = self.approx_max_mi(self.D, x, y)
            for i, v in enumerate(self.approx_max_mi(self.D, x, y)):
                I[i + 2][y] = v

            for i, v in enumerate(self.approx_max_mi(self.D_orth, x, y)):
                I_orth[i + 2][y] = v

        def characteristic_value(x, y):
            return max(I[x][y], I_orth[y][x]) if (x * y) <= b and x != 0 and y != 0 else np.nan

        I = np.fromfunction(np.vectorize(characteristic_value), (s, s), dtype=np.float64)

        def normalize(x, y):
            return I[x][y] / min(log2(x), log2(y)) if (x * y) <= b and x != 0 and y != 0 else np.nan

        M = np.fromfunction(np.vectorize(normalize), (s, s), dtype=np.float64)
        return M

    def create_partition(self, ordinals, axis='x'):
        assert axis == 'x' or axis == 'y'

        p = {}
        i = -1
        for p_start, p_end in pairwise(ordinals):
            i += 1
            for p_index in range(p_start + 1, p_end + 1):
                p[self.Dx[p_index]] = i
        return p

    def optimize_x_axis(self, d_x, q, x):
        c = self.get_c(self.get_clumps_partition(q))
        k = len(c) - 1
        I = np.zeros(shape=(k + 1, x + 1), dtype=np.float64)
        P = np.empty(shape=(k + 1, x + 1), dtype=np.ndarray)

        def F(s, t, l):
            return (c[s] / c[t]) * (I[s][l - 1] - self.hq(q)) - (c[t] - c[s] / c[t]) * self.hpq([c[s], c[t]], q)

        #Find the optimal partitions of size 2
        for t in xrange(2, k + 1):
            s = max(range(1, t + 1), key=lambda a: self.hp([c[0], c[a], c[t]]) - self.hpq([c[0], c[a], c[t]], q))
            P[t][2] = np.array([c[0], c[s], c[t]])
            I[t][2] = self.hq(q) + self.hp(P[t][2]) - self.hpq(P[t][2], q)

        #Inductively build the rest of the table of optimal partitions
        for l in xrange(3, x + 1):
            for t in xrange(l, k + 1):
                s = max(xrange(l - 1, t + 1), key=lambda s_: F(s_, t, l))
                P[t][l] = self.extend(P[s][l - 1], c[t])
                I[t][l] = self.hq(q) + self.hp(P[t][l]) - self.hpq(P[t][l], q)

        for l in range(k + 1, x + 1):
            I[k][l] = I[k][k]

        return I[k][2:x + 1]

    @staticmethod
    def equipartition_y_axis(d, y):
        n = len(d)

        desired_row_size = n / y
        i = 0
        sharp = 0
        current_row = 0
        q = {}
        while i < n:
            s = len([p for p in d if d[i][1] == p[1]])
            lhs = abs(float64(sharp) + float64(s) - desired_row_size)
            rhs = abs(float64(sharp) - desired_row_size)

            if sharp != 0 and lhs >= rhs:
                sharp = 0
                current_row += 1
                temp1 = n - i
                temp2 = y - current_row
                desired_row_size = temp1 / temp2

            for j in xrange(s):
                q[d[i + j]] = current_row

            i += s
            sharp += s

        return q

    def get_clumps_partition(self, q):
        q_tilde = copy(q)
        i = 0
        c = -1

        while i < self.n:
            s = 1
            flag = False
            for j in xrange(i + 1, self.n):
                if self.Dx[i][0] == self.Dx[j][0]:
                    s += 1
                    if q_tilde[self.Dx[i]] != q_tilde[self.Dx[j]]:
                        flag = True
                else:
                    break

            if s > 1 and flag:
                for j in xrange(s):
                    q_tilde[self.Dx[i + j]] = c
                c -= 1
            i += s

        i = 0
        p = {self.Dx[0]: 0}
        for j in xrange(1, self.n):
            if q_tilde[self.Dx[j]] != q_tilde[self.Dx[j - 1]]:
                i += 1
            p[self.Dx[j]] = i

        return p

    def get_super_clumps_partition(self, q, k_hat):
        p_tilde = self.get_clumps_partition(q)

        k = len(p_tilde)
        if k > k_hat:
            d_p_tilde = sorted([(0, p_tilde[p]) for p in self.Dx], key=lambda point: point[1])
            p_hat = self.equipartition_y_axis(d_p_tilde, k_hat)
            p = {point: p_hat[(0, p_tilde[point])] for point in self.Dx}
            return p
        else:
            return p_tilde

    def extend(self, ordinals, c):
        if any(c == existing for existing in ordinals):
            new_ordinals = ordinals
        else:
            from bisect import insort
            new_ordinals = list(ordinals)
            insort(new_ordinals, c)
        return np.array(new_ordinals)

    #TODO optimize
    def get_c(self, p_map):

        p_grouped = group_points_by_partition(p_map)
        c0 = [self.Dx.index(get_leftest_point(p_grouped[0])) - 1]

        def last_point_index(partition_index):
            return self.Dx.index(get_rightest_point(p_grouped[partition_index]))

        c1_k = map(last_point_index, p_grouped.iterkeys())
        c = c0 + c1_k
        return c

    def p_distr(self, ordinals):
        return np.fromiter((end - start for start, end in pairwise(ordinals)), dtype=int)

    def hp(self, ordinals):
        distribution = self.p_distr(ordinals)
        return entropy(distribution / distribution.sum())

    def hq(self, q):
        n = len(q)
        distribution = np.fromiter(Counter(q.itervalues()).itervalues(), dtype=int)
        return entropy(distribution / n)

    #TODO x_partition it would be better to be ordinals instead of map
    def hpq(self, x_ordinals, y_map):
        x_partition = self.create_partition(x_ordinals)
        grid_hist = self.get_grid_histogram(x_partition, y_map)
        return entropy(grid_hist / len(x_partition))

    def get_grid_histogram(self, p_map, q_map):
        rows, columns = group_points_by_partition(q_map), group_points_by_partition(p_map)

        def grid_cell_size(row_index, column_index):
            return len(set(rows[row_index]) & set(columns[column_index]))

        grid_points_distribution = (grid_cell_size(r, c) for r in reversed(xrange(len(rows))) for c in xrange(len(columns)))
        return np.fromiter(grid_points_distribution, dtype=int)


def entropy(P):
    '''
    Return the Shannon entropy of a probability vector P
    See http://www.scholarpedia.org/article/Entropy#Shannon_entropy
    '''
    h = -np.fromiter((i * np.log2(i) for i in P if i > 0), dtype=np.float64).sum()
    return h


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def group_points_by_partition(p):
    d = defaultdict(list)
    for k, v in p.iteritems():
        d[v].append(k)
    return dict(d)


def get_rightest_point(points):
    return max(points, key=p_x)


def get_leftest_point(points):
    return min(points, key=p_x)


def get_highest_point(points):
    return max(points, key=p_y)


def get_lowest_point(points):
    return min(points, key=p_y)


def last_abscissa(x_bin):
    return p_x(get_rightest_point(x_bin))


def last_ordinate(y_bin):
    return p_y(get_highest_point(y_bin))


def plot_partitions(p, q, file_name='example_grid.png', output_dir='/home/florents/workspace/mine/doc/examples/'):
    x_axis_partition, y_axis_partition = group_points_by_partition({point: p[point] for point in p.points()}), group_points_by_partition(q)

    from itertools import chain

    points = set(chain(p.points(), q.iterkeys()))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Scatter points
    ax.scatter(map(p_x, points), map(p_y, points))

    x_bin_edge = lambda x_bin: last_abscissa(x_bin) + 0.2
    y_bin_edge = lambda y_bin: last_ordinate(y_bin) + 0.2

    x_ticks = map(x_bin_edge, x_axis_partition.values())
    y_ticks = map(y_bin_edge, y_axis_partition.values())

    ax.get_xaxis().set_ticks(x_ticks)
    ax.get_yaxis().set_ticks(y_ticks)

    # Format grid appearance
    ax.grid(True, alpha=0.5, color='red', linestyle='--', linewidth=1.5)

    x_partition_size = len(x_axis_partition.values())
    y_partition_size = len(y_axis_partition.values())
    plt.title(str(x_partition_size) + ' - by - ' + str(y_partition_size) + ' Grid')
    plt.savefig(output_dir + file_name)


def plot_char_matrix_surface(m, file_name='example_grid.png', output_dir='/home/florents/workspace/mine/doc/examples/'):
    x = np.arange(m.shape[0])
    y = np.arange(m.shape[1])
    xx, yy = np.meshgrid(x, y)

    @np.vectorize
    def char_value(x, y):
        return m[x][y]

    z = char_value(xx, yy)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
    ax.set_zlim3d(-1.01, 1.01)

    fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.savefig(output_dir + file_name)
