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

"""
Each point is a tuple, thus:
P = (x, y) ==> 
P[0] = x and P[1] = y
"""

import bisect
from collections import defaultdict, Mapping
from copy import copy
from itertools import chain, imap
from math import log
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np


p_x, p_y = lambda p: p[0], lambda p: p[1]

def get_rightest_point(points):
    return max(points, key=lambda p: p[0])

def get_uppest_point(points):
    return max(points, key=lambda p: p[1])

def visualize_grid(x_axis_parition={}, y_axis_partition={}, step=0.2):
    """
    x_axis_parition: x-axis partition
    y_axis_partition: p-axis partition
    """
    points = set(chain(x_axis_parition.iterkeys(), y_axis_partition.iterkeys()))
    
    x_axis_parition = GroupPartitionsPoints(x_axis_parition)
    y_axis_partition = GroupPartitionsPoints(y_axis_partition)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.scatter([p[0] for p in points], [p[1] for p in points])
    
    x_ticks = [get_rightest_point(x_axis_parition[x_bin])[0] + step for x_bin in x_axis_parition.iterkeys()]
    y_ticks = [get_uppest_point(y_axis_partition[y_bin])[1] + step for y_bin in y_axis_partition.iterkeys()]
    if x_ticks: del x_ticks[-1]
    if y_ticks: del y_ticks[-1]
    
    
    ax.get_xaxis().set_ticks(x_ticks)
    ax.get_yaxis().set_ticks(y_ticks)
    ax.grid(True)
    plt.show()

def is_sorted_increasing_by(D, increasing_by='x'):
    assert increasing_by == 'x' or increasing_by == 'y'
    
    if increasing_by == 'x':
        return all(D[i][0] <= D[i + 1][0] for i in xrange(len(D) - 1))
    else:
        return all(D[i][1] <= D[i + 1][1] for i in xrange(len(D) - 1))

def sort_D_increasing_by(D, increasing_by='x'):
    assert increasing_by == 'x' or increasing_by == 'y'
    
    return sorted(D, key=p_x) if increasing_by == 'x' else sorted(D, key=p_y)
        
def visualize_partition_by_endpoint_indices(D, x_partition_endpoint_indices=[], y_partition_endpoint_indices=[], step=0.2):
    D_sorted_by_x = sort_D_increasing_by(D, 'x')
    D_sorted_by_y = sort_D_increasing_by(D, 'y')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter([p[0] for p in D], [p[1] for p in D])
    
    # Draw partition lines for each axis
    ax.get_xaxis().set_ticks([D_sorted_by_x[c][0] + step for c in x_partition_endpoint_indices])
    ax.get_yaxis().set_ticks([D_sorted_by_y[c][1] + step for c in y_partition_endpoint_indices])
   
    ax.grid(True)
    
    plt.show()
    
def visualize_partition_by_endpoints(D, x_endpoints, y_endpoints, step=0.2):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter([p[0] for p in D], [p[1] for p in D])
    
    # Draw partition lines for each axis
    ax.get_xaxis().set_ticks([p[0] + step for p in x_endpoints])
    ax.get_yaxis().set_ticks([p[1] + step for p in y_endpoints])
   
    ax.grid(True)
    
    plt.show()
    
def EquipartitionYAxis(D, y):
    if not is_sorted_increasing_by(D, 'y'): D = sort_D_increasing_by(D, 'y')
    
    n = len(D)
    
    desiredRowSize = float(n) / float(y)
    
    i = 0
    sharp = 0
    currRow = 0
    
    Q = {}
    while(i < n):
        #s = len([p for p in D if p[1] == D[i][1]])
        #s = sum(imap(lambda p: p[1] == D[i][1], D))
        
        s = sum(1 for p in D if p[1] == D[i][1])
        
        lhs = abs(float(sharp) + float(s) - desiredRowSize)
        rhs = abs(float(sharp) - desiredRowSize)
        
        if (sharp != 0 and lhs >= rhs):
            currRow += 1
            sharp = 0
            temp1 = float(n) - float(i)
            temp2 = float(y) - float(currRow)
            desiredRowSize = temp1 / temp2
        
        for j in range(s): 
            Q[D[i + j]] = currRow
        
        i += s
        sharp += s
    
    return Q

def GetClumpsPartition(D, Q):
    if not is_sorted_increasing_by(D, 'x'): D = sort_D_increasing_by(D, 'x')
    
    n = len(D)
    
    Q_tilde = copy(Q)
    
    i = 0
    c = -1 
    
    while(i < n):
        s = 1
        flag = False
        for j in range(i + 1, n):
            if D[i][0] == D[j][0]:
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
    if not is_sorted_increasing_by(D, 'x'): D = sort_D_increasing_by(D, 'x')

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
    if not is_sorted_increasing_by(D, 'x'): D = sort_D_increasing_by(D, 'x')
    
    super_clumps_partition = GetSuperclumpsPartition(D, Q, k_hat)
    c = GetPartitionOrdinals(D, super_clumps_partition, axis='x')
    #visualize_grid(super_clumps_partition)
    
    k = len(set(super_clumps_partition.values()))
    assert k == len(c) - 1

    #Find the optimal partitions of size 2 
    I = np.zeros(shape=(k+1, x+1))
    for t in range(1, len(c)):
        print t
        s = max(range(0,t), 
                key=lambda s: 
                            H(GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]], axis='x')) - 
                            H(GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]], axis='x'), Q))
        print [c[0], c[s], c[t]]
        P_t_2 = GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]], axis='x')
        print P_t_2
        I[t][2] = H(Q) + H(P_t_2) - H(P_t_2, Q)
    
    #Inductively build the rest of the table of optimal partitions
    for l in range(3, x+1):
        for t in range(l, k+1):
            s = max(range(l-1, t+1), 
                    key=lambda s: float((c[s]/c[t])) * (I[s][l-1] - H(Q)) - float(((c[t]-c[s]) / c[t])) * H(GetPartitionFromOrdinals(D, [c[s],c[t]], axis='x'), Q)
                    )
            ordinals_t_l_1 = [s, l-1]
            bisect.insort(ordinals_t_l_1, c[t])
            P_t_l = GetPartitionFromOrdinals(D, ordinals_t_l_1, axis='x')
            I[t][l] = H(Q) + H(P_t_l) - H(P_t_l, Q)
    
    #for l in range(k+1, x+1): 
    for l in range(k+1, x+1):
        I[k][l] = I[k,k]
    for l in range(k+1, x+1):
        I[k][l] = I[k][k]
    return [I[k][i] for i in range(2, x+1)]
            
        
def GetPartitionOrdinals(D, P, axis='x'):
    P_tilde = GroupPartitionsPoints(P)
    if axis == 'x':
        return [-1] + [D.index(get_rightest_point(P_tilde[k])) for k in sorted(P_tilde.keys())]
    elif axis == 'y':
        return [D.index(get_uppest_point(P_tilde[k])) for k in sorted(P_tilde.keys())]
  
  
from itertools import tee, izip
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def GetPartitionFromOrdinals(D, ordinals, axis='x'):
    if not is_sorted_increasing_by(D, 'x'): D = sort_D_increasing_by(D, 'x')
    P = {}
    
    current_partition = 0
    for i, j in pairwise(ordinals):
        from_point = i + 1
        to_point = j
        for p_index in range(from_point,to_point+1):
            P[D[p_index]] = current_partition
        
        current_partition += 1
    
    return P

def H(P, *Q):

    if not Q:
        return entropy(GetProbabilityDistribution(P)) if isinstance(P, Mapping) else entropy(P)
    
    else:
        # Not implemented for len(Q) > 1
        assert len(Q) == 1
        
        Q = Q[0]
        # Compute joint entropy
        n_points = float(len(set(chain(P.iterkeys(), Q.iterkeys()))))
        
        grid_matrix = GetGridMatrix(P, Q)
        probabilities = grid_matrix.flatten() / n_points
        
        return entropy(probabilities)

def entropy(probs): 
    return -sum(p * log(p, 2) for p in probs if p > 0)
                 
def I(P, Q):
    return H(P) + H(Q) + H(P,Q)

def GetPartitionEndpointIndices(partition, D, axis='x', step=0.3):
    assert axis == 'x' or axis == 'y'
    
    if axis == 'x' and not is_sorted_increasing_by(D, increasing_by='x'):
        D = sort_D_increasing_by(D, increasing_by='x')
    elif axis == 'y' and not is_sorted_increasing_by(D, increasing_by='y'):
        D = sort_D_increasing_by(D, increasing_by='y')
    
    d = GroupPartitionsPoints(partition)
    
    endpoint_indices = []
    for k in sorted(d.keys()):
        if axis == 'x':
            endpoint_x = max(d[k], key=lambda p: p[0])
            endpoint_indices.append(D.index(endpoint_x))
        elif axis == 'y':
            endpoint_y = max(d[k], key=lambda p: p[1])
            endpoint_indices.append(D.index(endpoint_y))

def GroupPartitionsPoints(P):
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

def GetProbabilityDistribution(P):
    """
    n: number of total points
    """
    n = len(set(P.keys()))
    d = GroupPartitionsPoints(P)    
    return [float(len(d[k])) / float(n) for k in sorted(d.keys())]

def GetGridMatrix(P, Q):
    """
    Each matrix element equals the number of points in the corresponding grid cell.
    """
    P = GroupPartitionsPoints(P)
    Q = GroupPartitionsPoints(Q)
   
    grid_matrix = np.zeros(shape=(len(Q.keys()), len(P.keys())), dtype=int)
    num_rows, num_columns = grid_matrix.shape[0], grid_matrix.shape[1]
    for r in range(num_rows):
        for c in range(num_columns):
            grid_matrix[r][c] = len(set(Q[r]) & set(P[c]))
    flipped = np.flipud(grid_matrix)
    return flipped
            
