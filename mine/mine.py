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


import bisect
from collections import defaultdict, Counter
from copy import copy
from itertools import *
from math import floor
from numpy import vstack, lexsort, shape, where, log2, fliplr

from computations import *
import matplotlib.pyplot as plt
import numpy as np
from utils import *


class MINE:
    def __init__(self,x,y):
        self.D = vstack((x,y)).T
        self.n = len(self.D)
        self.Dx_indices = lexsort((self.D[:,1],self.D[:,0]))
        self.Dx = self.D[self.Dx_indices]

        self.Dy_indices = lexsort((self.D[:,0],self.D[:,1]))
        self.Dy = self.D[self.Dy_indices]
        self.D_orth = fliplr(self.D)

    def get_point(self, index, axis_sorted_by):
        assert axis_sorted_by == 'x' or axis_sorted_by == 'y'

        return (self.Dx[index][0], self.Dx[index][1]) if axis_sorted_by=='x' else (self.Dy[index][0], self.Dy[index][1])

    def p_x(p):
        return p[0]

    def p_y(p):
        return p[1]

    @staticmethod
    def equipartition_y_axis(d_y, y):
        n = len(d_y)

        desired_row_size = float64(n) / float64(y)

        i = 0
        sharp = 0
        current_row = 0

        q = {}
        while i < n:
            s = shape(where(d_y[:,1] == d_y[i][1]))[1]
            lhs = abs(float64(sharp) + float64(s) - desired_row_size)
            rhs = abs(float64(sharp) - desired_row_size)

            if sharp != 0 and lhs >= rhs:
                sharp = 0
                current_row += 1
                temp1 = float64(n) - float64(i)
                temp2 = float64(y) - float64(current_row)
                desired_row_size = temp1 / temp2

            for j in xrange(s):
                point = (d_y[i+j][0], d_y[i+j][1])
                q[point] = current_row

            i += s
            sharp += s

        return q

    def get_points_assignments(self, d, axis_sorted_by='y'):
        return {self.get_point(k, axis_sorted_by): v for k,v in d.iteritems()}

    def get_clumps_partition(self, q):
        q_tilde = q.copy()
        i = 0
        c = -1

        while i < self.n:
            s = 1
            flag = False
            for j in xrange(i + 1, self.n):
                if p_x(self.get_point(i, 'x')) == p_x(self.get_point(j, 'x')):
                    s += 1
                    if q_tilde[self.get_point(i, 'x')] != q_tilde[self.get_point(j, 'x')]:
                        flag = True
                else:
                    break

            if s > 1 and flag:
                for j in xrange(s):
                    q_tilde[self.get_point(i+j, 'x')] = c
                c -= 1
            i += s

        i = 0
        p = {self.get_point(0, 'x'): 0}
        ordinals = [i-1]
        for j in xrange(1, self.n):
            if q_tilde[self.get_point(j, 'x')] != q_tilde[self.get_point(j-1, 'x')]:
                ordinals.append(j-1)
                i += 1
            if j == self.n-1:
                ordinals.append(j)
            p[self.get_point(j, 'x')] = i
        return p, ordinals

    def get_super_clumps_partition(self, q, k_hat):
        p_tilde, _ = self.get_clumps_partition(q)
        k = len(set(p_tilde.itervalues()))
        if k > k_hat:
            x = np.zeros((self.n,), dtype=np.int)
            y = np.fromiter(p_tilde.itervalues(), dtype=np.int)
            d_p_tilde = np.vstack((x,y)).T
            #Sort by increasing y-value
            d_p_tilde_y_indices = lexsort((d_p_tilde[:,0],d_p_tilde[:,1]))
            d_p_tilde = d_p_tilde[d_p_tilde_y_indices]
            p_hat = self.equipartition_y_axis(d_p_tilde, k_hat)
            p = {tuple(point):p_hat[(0, p_tilde[tuple(point)])] for point in self.D}
            print 'first'
            return p
        else:
            print 'in'
            return p_tilde


def get_all_size_2_partition(ordinals):
    k = len(ordinals)
    for t in xrange(2, k):
        for s in xrange(1, t+1):
            yield np.array((ordinals[0], ordinals[s], ordinals[t]), dtype=np.int32)


def HP(ordinals):
    #n = number_of_points_in_partition(ordinals)
    pass


def number_of_points_in_partition(ordinals):
    if ordinals[0] > 0:
        return sum(end_point - start_point for start_point, end_point in pairwise(ordinals))
    else:
        assert ordinals[0] == -1

        if len(ordinals) == 3:
            return 1 + ordinals[1] + (ordinals[2] - ordinals[1])

        return ordinals[1] - 1 + sum(end_point - start_point for start_point, end_point in pairwise(ordinals[1:]))


def get_partition_histogram(ordinals):
    return np.fromiter((end_point - start_point for start_point, end_point in pairwise(ordinals)), dtype=int)

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

def plot_partitions(p, q, file_name='example_grid.png', output_dir='/home/florents/workspace/mine/doc/examples/'):
    x_axis_partition, y_axis_partition = group_points_by_partition(p), group_points_by_partition(q)

    from itertools import chain
    points = set(chain(p.iterkeys(), q.iterkeys()))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Scatter points
    ax.scatter(map(p_x, points), map(p_y, points))

    x_bin_edge = lambda x_bin :last_abscissa(x_bin) + 0.2
    y_bin_edge = lambda y_bin:  last_ordinate(y_bin) + 0.2

    x_ticks = map(x_bin_edge, x_axis_partition.values())
    y_ticks = map(y_bin_edge, y_axis_partition.values())

    ax.get_xaxis().set_ticks(x_ticks)
    ax.get_yaxis().set_ticks(y_ticks)

    # Format grid appearance
    ax.grid(True, alpha=0.5, color='red', linestyle='--', linewidth=1.5)

    x_partition_size = len(x_axis_partition.values())
    y_partition_size = len(y_axis_partition.values())
    plt.title(str(x_partition_size) + ' - by - ' + str(y_partition_size) + ' Grid')
    plt.savefig(output_dir+file_name)

def ApproxMaxMI(D, x, y, k_hat):
    assert x > 1 and y > 1 and k_hat > 1 
    
    Q = EquipartitionYAxis(sort_D_increasing_by(D, increasing_by='y'), y)
    D = sort_D_increasing_by(D, increasing_by='x')
    return OptimizeXAxis(D, Q, x, k_hat)

def ApproxCharacteristicMatrix(D, B, c):
    assert B > 3 and c > 0
    
    D_orth = [tuple(reversed(p)) for p in D]
    
    s = int(floor(B / 2)) + 1
    
    I = np.zeros(shape=(s, s), dtype=float64)
    I_orth = np.zeros(shape=(s, s), dtype=float64)
    M = np.zeros(shape=(s, s), dtype=float64)

    '''
    Lines 2-6
    '''
    for y in xrange(2, s):
        x = int(floor(B / y))
        
        for i, v in enumerate(ApproxMaxMI(D, x, y, c * x)): I[i + 2][y] = v
 
        for i, v in enumerate(ApproxMaxMI(D_orth, x, y, c * x)): I_orth[i + 2][y] = v
        
        
    '''
    Lines 7-10
    '''
    def characteristic_value(x, y):
        return max(I[x][y], I_orth[y][x]) if (x * y) <= B and x != 0 and y != 0 else np.nan
        
    I = np.fromfunction(np.vectorize(characteristic_value), (s, s), dtype=np.float64)
    
    def normalize(x, y):
        return I[x][y] / min(log2(x), log2(y)) if (x * y) <= B and x != 0 and y != 0 else np.nan
        
    M = np.fromfunction(np.vectorize(normalize), (s, s), dtype=np.float64)
   
    return M
    
def OptimizeXAxis(D, Q, x, k_hat):
    assert is_sorted_increasing_by(D, 'x')
    
    super_clumps_partition = GetSuperclumpsPartition(D, Q, k_hat)
    c = GetPartitionOrdinalsFromMap(D, super_clumps_partition, axis='x')
    
    # Total number of clumps
    k = len(set(super_clumps_partition.itervalues()))
    assert k == len(c) - 1

    # Find the optimal partitions of size 2 
    I = np.zeros(shape=(k + 1, x + 1))
    
    for t in xrange(2, k + 1):
        s = max(xrange(1, t + 1),
                key=lambda s: HP([c[0], c[s], c[t]]) - HPQ([c[0], c[s], c[t]], Q))
        
        # Optimal partition of size 2 on the first t clumps
        P_t_2 = [c[0], c[s], c[t]]
        I[t][2] = HQ(Q) + HP(P_t_2) - HPQ(P_t_2, Q)
   
    # Inductively build the rest of the table of optimal partitions
    for l in xrange(3, x + 1):
        for t in xrange(l, k + 1):
            s = max(xrange(l - 1, t + 1),
                    key=lambda s: float64((c[s] / c[t])) * (I[s][l - 1] - HQ(Q)) - float64(((c[t] - c[s]) / c[t])) * HPQ([c[s], c[t]], Q)
                    )
            P_t_l = c[1:l - 1]
            bisect.insort(P_t_l, c[t])
            # Optimal partition of size l on the first t clumps of D
            I[t][l] = HQ(Q) + HP(P_t_l) - HPQ(P_t_l, Q)

    for l in xrange(k + 1, x + 1): I[k][l] = I[k][k]
    
    return I[k][2:x + 1]


def mine(cm, B, e=1):
    print np.indices(cm.shape)
    mic = max(value for (x, y), value in np.ndenumerate(cm) if x * y < B and (x, y) != (0, 0) and not np.isnan(value))
    mas = max(np.abs(value - cm[y][x]) for (x, y), value in np.ndenumerate(cm) if x * y < B and (x, y) != (0, 0) and not np.isnan(value))
    mev = max(value for (x, y), value in np.ndenumerate(cm) if x * y < B and (x, y) != (0, 0) and not np.isnan(value) and (x is 2 or y is 2))
    mcn = min(np.log2(x * y) for (x, y), value in np.ndenumerate(cm) if x * y < B and (x, y) != (0, 0) and not np.isnan(value) and value >= (1 - e) * mic)
    
    return {'MIC':mic,
            'MAS':mas,
            'MEV':mev,
            'MCN':mcn
            }
def mic(x,y):
    D = zip(x, y)
    n = len(D)
    B = pow(n, 0.6)
    c = 15
    M = ApproxCharacteristicMatrix(D, B, c=1)
    return mine(M, B, c)['MIC']