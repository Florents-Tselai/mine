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
from collections import defaultdict, Counter
from copy import copy
from itertools import *
from math import floor

from numpy import vstack, lexsort, shape, where, log2, fliplr, float64
import matplotlib.pyplot as plt
import numpy as np


p_x, p_y = lambda p: p[0], lambda p: p[1]


class MINE:
    def __init__(self, x, y):
        self.D = np.core.records.fromarrays([x,y], names='x,y')
        self.D_orth = np.core.records.fromarrays([y,x], names='x,y')
        self.Dx = self.D[self.D.argsort(order='x')]
        self.Dy = self.D[self.D.argsort(order='y')]
        self.n = len(self.D)

    def create_partition(self, ordinals, axis='x'):
        assert axis == 'x' or axis=='y'

        d = self.Dx if axis=='x' else self.Dy
        return Partition(d=d, ordinals=ordinals)

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

        for y in range(2, int(floor(b / 2))+1):
            x = int(floor(b / y))
            assert len(self.approx_max_mi(self.D, x, y)) == x-1
            #print len(I[2:x+1][y]), len(self.approx_max_mi(self.D, x, y))
            #I[2:x+1][y] = self.approx_max_mi(self.D, x, y)
            #for i, v in enumerate(self.approx_max_mi(self.D, x, y)): I[i + 2][y] = v

            for i, v in enumerate(self.approx_max_mi(self.D_orth, x, y)): I_orth[i + 2][y] = v

        def characteristic_value(x, y):
            return max(I[x][y], I_orth[y][x]) if (x * y) <= b and x != 0 and y != 0 else np.nan

        I = np.fromfunction(np.vectorize(characteristic_value), (s, s), dtype=np.float64)

        def normalize(x, y):
            return I[x][y] / min(log2(x), log2(y)) if (x * y) <= b and x != 0 and y != 0 else np.nan

        M = np.fromfunction(np.vectorize(normalize), (s, s), dtype=np.float64)
        return M

    def optimize_x_axis(self, d_x, q, x):
        c = self.get_clumps_partition(q).ordinals
        k = len(c) - 1
        I = np.zeros(shape=(k + 1, x + 1), dtype=np.float64)
        P = np.empty(shape=(k+1, x+1), dtype=object)

        for t in xrange(2, k + 1):
            s = max(range(1, t + 1), key=lambda a: self.hp([c[0], c[a], c[t]]) - self.hpq(self.create_partition([c[0], c[a], c[t]]), q))
            P[t][2] = self.create_partition([c[0], c[s], c[t]])
            I[t][2] = q.h() + P[t][2].h() - self.hpq(P[t][2], q)
        for l in xrange(3, x + 1):
            for t in xrange(l, k + 1):
                def f(s_, t_, l_):
                    return (c[s_] / c[t_]) * (I[s_][l_ - 1] - q.h()) - (
                        ((c[t_] - c[s_] / c[t_])) * self.hpq(self.create_partition([c[s_], c[t_]]), q))

                s = max(xrange(l - 1, t + 1), key=lambda a: f(a, t, l))
                P[t][l] = P[t][l-1]+c[t]
                I[t][l] = q.h() + P[t][l].h() - self.hpq(P[t][l], q)
        for l in range(k + 1, x + 1):
            I[k][l] = I[k][k]

        return I[k][2:x + 1]

    @staticmethod
    def equipartition_y_axis(d_y, y):
        n = len(d_y)

        desired_row_size = float64(n) / float64(y)

        i = 0
        sharp = 0
        current_row = 0
        q = {}
        while i < n:
            s = shape(where(d_y['y'] == d_y[i]['y']))[1]
            lhs = abs(float64(sharp) + float64(s) - desired_row_size)
            rhs = abs(float64(sharp) - desired_row_size)

            if sharp != 0 and lhs >= rhs:
                sharp = 0
                current_row += 1
                temp1 = float64(n) - float64(i)
                temp2 = float64(y) - float64(current_row)
                desired_row_size = temp1 / temp2

            for j in xrange(s):
                point = (d_y[i + j][0], d_y[i + j][1])
                q[point] = current_row

            i += s
            sharp += s

        return Partition(map_assignments=q)

    def get_clumps_partition(self, q):
        q_tilde = copy(q)
        i = 0
        c = -1

        while i < self.n:
            s = 1
            flag = False
            for j in xrange(i + 1, self.n):
                if self.Dx[i]['x'] == self.Dx[j]['x']:
                    s += 1
                    if q_tilde[tuple(self.Dx[i])] != q_tilde[tuple(self.Dx[j])]:
                        flag = True
                else:
                    break

            if s > 1 and flag:
                for j in xrange(s):
                    q_tilde[tuple(self.Dx[i+j])] = c
                c -= 1
            i += s

        i = 0
        p = {tuple(self.Dx[0]): 0}
        ordinals = [i - 1]
        for j in xrange(1, self.n):
            if q_tilde[tuple(self.Dx[j])] != q_tilde[tuple(self.Dx[j-1])]:
                ordinals.append(j - 1)
                i += 1
            if j == self.n - 1:
                ordinals.append(j)
            p[tuple(self.Dx[j])] = i

        return Partition(d=self.Dx, ordinals=ordinals)

    def get_super_clumps_partition(self, q, k_hat):
        p_tilde= self.get_clumps_partition(q)
        k = len(p_tilde)
        if k > k_hat:
            x = np.arange(self.n)
            y = np.fromiter((p_tilde[tuple(p)] for p in self.Dx), dtype=np.int32)
            d_p_tilde = np.core.records.fromarrays([x,y], names='x,y')
            p_hat = self.equipartition_y_axis(d_p_tilde, k_hat)
            b = [set() for i in range(len(set(p_hat.values())))]

            for point_index, bin_index in p_hat.iteritems():
                b[bin_index].add(tuple(self.Dx[point_index[0]]))
            part = Partition(None, None, b)
            return part

        else:
            return p_tilde

    def hp(self, ordinals):
        return self.create_partition(ordinals=ordinals).h()

    def hq(self, q):
        return q.h()

    def hpq(self, x_partition, y_map):
        grid_hist = x_partition.grid_histogram(y_map)
        return entropy(grid_hist / x_partition.number_of_points())


class Partition:
    def __init__(self, d=None, ordinals=None,map_assignments=None):
        self.ordinals = ordinals
        if map_assignments is None:
            self.d = d
            self.bins = [set(self._get_point(i) for i in xrange(start+1, end+1)) for start, end in pairwise(ordinals)]
            self.map_assignments = {}
            for i, b in enumerate(self.bins):
                for p in b:
                    self.map_assignments[p] = i
        else:
            self.map_assignments = map_assignments
            b = [set() for i in range(len(set(map_assignments.values())))]

            for point, bin_index in map_assignments.iteritems():
                    b[bin_index].add(point)
            self.bins = b

    def _get_point(self, i):
        point = (self.d[i][0], self.d[i][1])
        return point

    def __len__(self):
        return len(self.bins)

    def number_of_points(self):
        return self.histogram().sum()

    def histogram(self):
        return np.fromiter((len(b) for b in self),dtype=np.int32)

    def grid_histogram(self, q):
        rows, columns = group_points_by_partition(q.map_assignments), self.bins

        def grid_cell_size(row_index, column_index):
            return len(set(rows[row_index]) & set(columns[column_index]))

        grid_points_distribution = (grid_cell_size(r, c) for r in reversed(xrange(len(rows))) for c in xrange(len(columns)))
        return np.fromiter(grid_points_distribution, dtype=int)


    def __copy__(self):
        return copy(self.map_assignments)

    def __getitem__(self, p):
        return self.map_assignments[p]

    def __iter__(self):
        return iter(self.bins)

    def h(self):
        return entropy(self.histogram() / self.number_of_points())

    def points(self):
        return frozenset().union(*self.bins)

    def __add__(self, c):
        assert self.ordinals is not None
        if any(c==existing for existing in self.ordinals):
            new_ordinals = self.ordinals
        else:
            from bisect import insort
            new_ordinals = list(self.ordinals)
            insort(new_ordinals, c)
        p = Partition(d=self.d, ordinals=new_ordinals)
        return p

    def __str__(self):
        return str(self.map_assignments)

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
    x_axis_partition, y_axis_partition = group_points_by_partition({point:p[point] for point in p.points()}), group_points_by_partition(q)

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

