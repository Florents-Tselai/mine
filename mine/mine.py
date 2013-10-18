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

import bisect
from collections import defaultdict, Mapping
from copy import copy
from itertools import chain, tee, izip
from math import log, floor

import matplotlib.pyplot as plt
import numpy as np


p_x, p_y = lambda p: p[0], lambda p: p[1]

def get_rightest_point(points): return max(points, key=p_x)

def get_leftest_point(points): return min(points, key=p_x)

def get_uppest_point(points): return max(points, key=p_y)

def get_downest_point(points): return min(points, key=p_y)

def last_abscissa(x_bin): return p_x(get_rightest_point(x_bin))
def last_ordinate(y_bin): return p_y(get_uppest_point(y_bin))

def ApproxMaxMI(D, x, y, k_hat):
    assert x > 1 and y > 1 and k_hat > 1 
    
    Q = EquipartitionYAxis(sort_D_increasing_by(D, increasing_by='y'), y)
    D = sort_D_increasing_by(D, increasing_by='x')
    return OptimizeXAxis(D, Q, x, k_hat)

def ApproxCharacteristicMatrix(D, B, c):
    assert B > 3 and c > 0
    
    D_orth = [tuple(reversed(p)) for p in D]
    
    s = int(floor(B / 2)) + 1
    
    I = np.zeros(shape=(s, s))
    I_orth = np.zeros(shape=(s, s))
    M = np.zeros(shape=(s, s))

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
    for x in xrange(2, s):
        for y in xrange(2, s):
            if x * y > B:
                continue
            else:
                I[x][y] = max(I[x][y], I_orth[y][x])
                M[x][y] = float(I[x][y]) / min(log(x), log(y))
    return M
    
def EquipartitionYAxis(D, y):
    assert is_sorted_increasing_by(D, 'y')
    
    n = len(D)
    
    desiredRowSize = float(n) / float(y)
    
    i = 0
    sharp = 0
    currRow = 0
    
    Q = {}
    while(i < n):
        # s = len([p for p in D if p[1] == D[i][1]])
        # s = sum(imap(lambda p: p[1] == D[i][1], D))
        
        s = sum(1 for p in D if p[1] == D[i][1])
        
        lhs = abs(float(sharp) + float(s) - desiredRowSize)
        rhs = abs(float(sharp) - desiredRowSize)
        
        if (sharp != 0 and lhs >= rhs):
            sharp = 0
            currRow += 1
            temp1 = float(n) - float(i)
            temp2 = float(y) - float(currRow)
            desiredRowSize = temp1 / temp2
        
        for j in range(s): 
            Q[D[i + j]] = currRow
        
        i += s
        sharp += s
    
    return Q

def GetClumpsPartition(D, Q):
    assert is_sorted_increasing_by(D, 'x')
    
    n = len(D)
    
    Q_tilde = copy(Q)
    
    i = 0
    c = -1 
    
    while(i < n):
        s = 1
        flag = False
        for j in range(i + 1, n):
            if p_x(D[i]) == p_x(D[j]):
                s += 1
                if Q_tilde[D[i]] != Q_tilde[D[j]]:
                    flag = True
            else:
                break
            
        if s > 1 and flag:
            for j in range(s):
                Q_tilde[D[i + j]] = c
            c -= 1
        i += s
    
    i = 0
    P = {}
    P[D[0]] = 0
    for j in range(1, n):
        if Q_tilde[D[j]] != Q_tilde[D[j - 1]]:
            i += 1
        P[D[j]] = i
    
    return P

def GetSuperclumpsPartition(D, Q, k_hat):
    assert is_sorted_increasing_by(D, 'x')

    P_tilde = GetClumpsPartition(D, Q)
    k = len(set(P_tilde.values()))    
    if k > k_hat:
        D_P_tilde = [(0, P_tilde[p]) for p in D]
        P_hat = EquipartitionYAxis(D_P_tilde, k_hat)
        P = {p:P_hat[(0, P_tilde[p])] for p in D}
        return P
    else:
        return P_tilde    
    
def OptimizeXAxis(D, Q, x, k_hat):
    assert is_sorted_increasing_by(D, 'x')
    
    super_clumps_partition = GetSuperclumpsPartition(D, Q, k_hat)
    c = GetPartitionOrdinalsFromMap(D, super_clumps_partition, axis='x')
    
    # Total number of clumps
    k = len(set(super_clumps_partition.values()))
    assert k == len(c) - 1

    # Find the optimal partitions of size 2 
    I = np.zeros(shape=(k + 1, x + 1))
    
    for t in xrange(2, k + 1):
        s = max(xrange(1, t + 1),
                key=lambda s: HP(D, [c[0], c[s], c[t]]) - HPQ([c[0], c[s], c[t]], Q))
        
        # Optimal partition of size 2 on the first t clumps
        P_t_2 = [c[0], c[s], c[t]]
        I[t][2] = HQ(Q) + HP(D, P_t_2) - HPQ(P_t_2, Q)
   
    # Inductively build the rest of the table of optimal partitions
    for l in xrange(3, x + 1):
        for t in xrange(l, k + 1):
            s = max(range(l - 1, t + 1),
                    key=lambda s: float((c[s] / c[t])) * (I[s][l - 1] - HQ(Q)) - float(((c[t] - c[s]) / c[t])) * HPQ([c[s], c[t]], Q)
                    )
            P_t_l = c[1:l - 1]
            bisect.insort(P_t_l, c[t])
            # Optimal partition of size l on the first t clumps of D
            I[t][l] = HQ(Q) + HP(D, P_t_l) - HPQ(P_t_l, Q)

    for l in xrange(k + 1, x + 1): I[k][l] = I[k][k]
    
    return I[k][2:x + 1]

def entropy(probs): 
    return -sum(p * log(p, 2) for p in probs if p > 0)
                 
def I(P_ordinals, Q_map):
    return HQ(Q_map) + HP(P_ordinals) - HPQ(P_ordinals, Q_map)

def HP(Dx, P_ordinals):
    assert is_sorted_increasing_by(Dx, 'x')
    
    # Number of points in the partition 
    n = m(P_ordinals)
    return entropy(np.array(GetPartitionHistogram(Dx, P_ordinals)) / float(n))

def HQ(Q_map):
    n = len(Q_map)
    temp = GroupPointsByPartition(Q_map)
    return entropy([len(temp[k]) / float(n) for k in temp])

def HPQ(P_ordinals, Q_map):
    Dx = sort_D_increasing_by(Q_map.keys(), 'x')
    
    return entropy(np.array(GetGridHistogram(Q_map, P_ordinals)) / float(m(P_ordinals)))

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def is_sorted_increasing_by(D, increasing_by='x'):
    assert increasing_by == 'x' or increasing_by == 'y'
    
    if increasing_by == 'x':
        return all(p_x(D[i]) <= p_x(D[i + 1]) for i in xrange(len(D) - 1))
    else:
        return all(p_y(D[i]) <= p_y(D[i + 1]) for i in xrange(len(D) - 1))

def m(ordinals):
    return sum(o2+1 if o1<0 else o2-o1 for o1, o2 in pairwise(ordinals))

def sort_D_increasing_by(D, increasing_by='x'):
    assert increasing_by == 'x' or increasing_by == 'y'
    
    return sorted(D, key=p_x) if increasing_by == 'x' else sorted(D, key=p_y)
       

def GetPartitionOrdinalsFromMap(D, P, axis='x'):
    assert is_sorted_increasing_by(D, axis)
    
    P_tilde = GroupPointsByPartition(P)
    
    if axis == 'x':
        return [D.index(get_leftest_point(P_tilde[0])) - 1] + [D.index(get_rightest_point(P_tilde[k])) for k in sorted(P_tilde.keys())]
    elif axis == 'y':
        return [D.index(get_downest_point(P_tilde[0])) - 1] + [D.index(get_uppest_point(P_tilde[k])) for k in sorted(P_tilde.keys())]
  
def GetPartitionMapFromOrdinals(D, ordinals, axis='x'):
    assert is_sorted_increasing_by(D, axis)
    
    to_be_binned = range(len(D))
        
    # Translate Reshef's convention to adhere to Numpy's one
    bins = [o + 1 for o in ordinals[:-1]]  
    # Assign point indices to bins formed by the partition ordinals
    map = {D[point_index]:partition - 1 for point_index, partition in enumerate(np.digitize(to_be_binned, bins))}
    
    return map

def GroupPointsByPartition(P):
    """
    P : point -> Partition index
    Returns
    d : partition index -> points
    
    Example:
    P = 
    {
    p1 -> 1
    p2 -> 2
    P3 -> 1
    p4 -> 2
    }
    
    Returns
    d = 
    {
    1 -> [p1, p3]
    2 -> [p2, p4]
    }
      
    """
    d = defaultdict(list)
    for k, v in P.iteritems(): 
        d[v].append(k)
    return dict(d)

def GetPartitionHistogram(D, ordinals, axis='x'):
    assert is_sorted_increasing_by(D, axis)
    return [p[1] + 1 if p[0] < 0 else p[1] - p[0] for p in pairwise(ordinals)]

def visualize(x_axis_parition={}, y_axis_partition={}, step=0.2):
    points = set(chain(x_axis_parition.iterkeys(), y_axis_partition.iterkeys()))
    
    x_axis_parition = GroupPointsByPartition(x_axis_parition)
    y_axis_partition = GroupPointsByPartition(y_axis_partition)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # Scatter points
    ax.scatter(map(p_x, points), map(p_y, points))
    
    x_bin_edge = lambda x_bin :last_abscissa(x_bin) + step
    y_bin_edge = lambda y_bin:  last_ordinate(y_bin) + step
    
    x_ticks = map(x_bin_edge, x_axis_parition.itervalues())
    y_ticks = map(y_bin_edge, y_axis_partition.itervalues())
    
    ax.get_xaxis().set_ticks(x_ticks)
    ax.get_yaxis().set_ticks(y_ticks)
    
    # Format grid appearance
    ax.grid(True, alpha=0.5, color='red', linestyle='-', linewidth=1.5)
    
    x_partition_size = len(x_axis_parition.values())
    y_partition_size = len(y_axis_partition.values())
    plt.title(str(x_partition_size) + ' - by - ' + str(y_partition_size) + ' Grid')
    plt.show()

def GetGridHistogram(Q, P_ordinals):
    rows = GroupPointsByPartition(Q)
    
    Dx = sort_D_increasing_by(Q.keys(), 'x')
    
    def cell_size(p1, p2, r): return sum(1 for point in Dx[(p1 + 1):(p2 + 1)] if point in rows[r])
    
    return [cell_size(p1, p2, r) for p1, p2 in pairwise(P_ordinals) for r in xrange(len(rows))]
