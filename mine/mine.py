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
import sys

p_x, p_y = lambda p: p[0], lambda p: p[1]


class MINE:
    def __init__(self, x, y, alpha=0.6, c=15):
        self.D = np.core.records.fromarrays([x, y], names='x,y')
        self.D_orth = list(np.core.records.fromarrays([y, x], names='x,y'))
        self.Dx = [tuple(p) for p in self.D[self.D.argsort(order='x')]]
        self.Dy = [tuple(p) for p in self.D[self.D.argsort(order='y')]]
        self.n = len(self.D)
        self.alpha = alpha
        self.c = c;
        self.b = int(self.n ** self.alpha)

    def compute(self):
        self.char_matrix = self.approx_char_matrix(self.D, self.b, self.c)

    def mic(self):
        return self.char_matrix.max()

    def approx_max_mi(self, d, x, y, k_hat):
        q = self.equipartition_y_axis(self.Dy, y)
        return self.optimize_x_axis(self.Dx, q, x, k_hat)

    def approx_char_matrix(self, d, b, c):
        from math import ceil, floor

        s = int(floor(b/2))

        I = np.zeros(shape=(s+1, s+1), dtype=float64)
        I_orth = np.zeros(shape=(s+1, s+1), dtype=float64)
        M = np.zeros(shape=(s+1, s+1), dtype=float64)

        for y in xrange(2, int(floor(b / 2)) + 1):
            x = int(floor(b / y))

            for i, v in enumerate(self.approx_max_mi(self.D, x, y, c*x)):
                I[i + 2][y] = v

            for i, v in enumerate(self.approx_max_mi(self.D_orth, x, y, c*x)):
                I_orth[i + 2][y] = v

        for x in xrange(2, s+1):
            for y in xrange(2, s+1):
                if x*y > b:
                    continue
                I[x][y] = max(I[x][y], I_orth[y][x])
                M[x][y] = I[x][y] / min(log2(x), log2(y))

        M = np.nan_to_num(M)
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

    def compute_cumhist(self, q_map, p_map):
        q = len(set(q_map.itervalues()))
        p = len(set(p_map.itervalues()))
        print q, p
        cumhist = np.zeros(shape=(q, p), dtype=int)

        for i in range(self.n):
            cumhist[q_map[self.Dx[i]]][p_map[self.Dx[i]]] += 1

        for i in range(q):
            for j in range(1, p):
                cumhist[i][j] += cumhist[i][j-1]
        return cumhist

    def hq(self, cumhist, q, p, n):
        prob, H = 0.0
        total = float(n)

        for i in range(q):
            if cumhist[i][p-1] != 0:
                prob = cumhist[i][p-1] / total
                H -= prob * log2(prob)
        return H

    def hp3(self, c, s, t):
        sum = 0
        prob, H = 0.0

        total = float(c[t-1])

        if s==t:
            return 0.0
        if c[s-1] != 0:
            prob = c[s-1] / total
            H -= prob*log2(prob)
        sum = c[t-1] - c[s-1]

        if sum != 0:
            prob = sum / total
            H -= prob *log2(prob)
        return H

    def hp3q(self, cumhist, c, q, p, s, t):
        prob, H = 0.0
        total = float(c[t-1])

        for i in range(q):
            if cumhist[i][s-1] != 0:
                prob = cumhist[i][s-1] / total
                H -= prob * log2(prob)
            sum = cumhist[i][t-1] - cumhist[i][s-1]
            if sum != 0:
                prob = sum / total
                H -= prob * log2(prob)
        return H

    def hp2q(self, cumhist, c, q, p, s, t):
        prob, H = 0.0, 0.0
        total = float(c[t-1] - c[s-1])

        if s == t:
            return 0.0
        for i in range(q):
            sum = cumhist[i][t-1] - cumhist[i][s-1]
            if sum != 0:
                prob = sum / total
                H -= prob * log2(prob)
        return H

    def compute_HP2Q(self, cumhist, c, q, p):
        HP2Q = np.zeros(shape=(p+1, p+1), dtype=int)
        for t in range(3, p+1):
            for s in range(2, t+1):
                HP2Q[s][t] = self.hp2q(cumhist, c, q, p, s, t)


    def optimize_x_axis(self, d_x, q, x, k_hat=sys.maxint):
        p_map = self.get_super_clumps_partition(q, k_hat)
        q_map = q
        p = len(set(p_map.iterkeys()))
        q = len(set(q_map.iterkeys()))

        if p == 1:
            return np.array([0 for i in range(x-1)], dtype=float)

        c = self.compute_c(p_map)
        #print c
        cumhist = self.compute_cumhist(q_map, p_map)
        HP2Q = self.compute_HP2Q(cumhist, c, q, p)
        HQ = self.hq(cumhist, q, p, n)
        print cumhist
        k = len(c) - 1
        I = np.zeros(shape=(k + 1, x + 1), dtype=np.float64)
        P = np.empty(shape=(k + 1, x + 1), dtype=np.ndarray)

        def F(s, t, l):
            return (c[s] / c[t]) * (I[s][l - 1] - self.hq(q)) - (c[t] - c[s] / c[t]) * self.hpq([c[s], c[t]], q)

        #Find the optimal partitions of size 2
        for t in xrange(2, p + 1):
            F_max = - np.inf
            for s in range(1, t+1):
                F = self.hp3(c, s, t) - self.hp3q(cumhist, c, q, p, s, t)
                if F > F_max:
                    I[t][2] = HQ + F
                    F_max = F

        #Inductively build the rest of the table of optimal partitions
        for l in xrange(3, x + 1):
            for t in xrange(1, p + 1):
                ct = float(c[t-1])
                F_max = -np.inf
                for s in range(l-1, t+1):
                    cs = float(c[s-1])
                    F = (cs/ct * I[s][l-1]-HQ) - (((ct-cs)/ct) * HP2Q[s][t])
                    if F > F_max:
                        I[t][l] = HQ + F
                        F_max = F

        for i in range(p + 1, x + 1):
            I[p][i] = I[p][p]

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
            d_p_tilde = [(0, p_tilde[p]) for p in self.Dx]
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
    def compute_c(self, p_map):
        p = len(set(p_map.itervalues()))
        c = [0 for i in range(p)]
        for i in range(self.n):
            c[p_map[self.Dx[i]]] += 1
        for i in range(1, p):
            c[i] += c[i-1]
        assert len(c) == p
        return c

    def p_distr(self, ordinals):
        return np.fromiter((end - start for start, end in pairwise(ordinals)), dtype=int)

    def hp(self, ordinals):
        distribution = self.p_distr(ordinals)
        return entropy(distribution / distribution.sum())

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
