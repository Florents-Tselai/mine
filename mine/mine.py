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
from itertools import chain, tee, izip
from math import log, floor
import matplotlib.pyplot as plt
import numpy as np

'''
Algorithms
'''
'''
Algorithm 4
'''
def ApproxMaxMI(D, x, y, k_hat):
    assert x>1 and y>1 and k_hat >1 
    
    Q = EquipartitionYAxis(D, y)
    D = sort_D_increasing_by(D, increasing_by='x')
    return OptimizeXAxis(D, Q, x, k_hat)

'''
Algorithm 5
'''
def ApproxCharacteristicMatrix(D, B, c):
    assert B > 3 and c > 0
    
    D_orth = [tuple(reversed(p)) for p in D]
    
    I = np.zeros(shape=(2,int(floor(B/2))))
    I_orth = np.zeros(shape=(2,int(floor(B/2))))
    M = np.zeros(shape=(int(floor(B/2))+1,int(floor(B/2))+1))
    
    '''
    Lines 2-6
    '''
    for y in range(2, int(floor(B/2))+1):
        x = int(floor(B/y))
        temp1 = ApproxMaxMI(D, x, y, c*x)
        
        I = np.append(I, temp1, axis=0)
        temp2 = ApproxMaxMI(D_orth, x, y, c*x)
        I_orth = np.append(I_orth, temp2, axis=0)
        
    '''
    Lines 7-10
    '''
    for x in range(2,int(floor(B/2))+1):
        for y in range(2,int(floor(B/2))+1):
            if x*y > B:
                continue
            else:
                I[x][y] = max(I[x][y], I_orth[y][x])
                M[x][y] = float(I[x][y]) / min(log(x), log(y))
    return M
    
'''
Algorithm 3
'''
def EquipartitionYAxis(D, y):
    if not is_sorted_increasing_by(D, 'y'): D = sort_D_increasing_by(D, 'y')
    
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
    if not is_sorted_increasing_by(D, 'x'): D = sort_D_increasing_by(D, 'x')
    
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
    
'''
Algorithm 2
'''
def OptimizeXAxis(D, Q, x, k_hat):
    if not is_sorted_increasing_by(D, 'x'): D = sort_D_increasing_by(D, 'x')
    
    super_clumps_partition = GetSuperclumpsPartition(D, Q, k_hat)
    c = GetPartitionOrdinals(D, super_clumps_partition, axis='x')
    
    # Total number of clumps
    k = len(set(super_clumps_partition.values()))
    assert k == len(c) - 1

    # Find the optimal partitions of size 2 
    I = np.zeros(shape=(k + 1, x + 1))
    
    for t in range(2, k+1):
        s = max(range(1, t+1),
                key=lambda s: 
                            H(GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]], axis='x')) - 
                            H(GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]], axis='x'), Q))
        # Optimal partition of size 2 on the first t clumps
        optimal_2_partition_ordinals =  [c[0], c[s], c[t]]
        P_t_2 = GetPartitionFromOrdinals(D, ordinals=optimal_2_partition_ordinals, axis='x')
        I[t][2] = H(Q) + H(P_t_2) - H(P_t_2, Q)
   
    # Inductively build the rest of the table of optimal partitions
    for l in range(3, x + 1):
        for t in range(l, k + 1):
            
            s = max(range(l - 1, t + 1),
                    key=lambda s: float((c[s] / c[t])) * (I[s][l - 1] - H(Q)) - float(((c[t] - c[s]) / c[t])) * H(GetPartitionFromOrdinals(D, [c[s], c[t]], axis='x'), Q)
                    )
            ordinals_t_l_1 = [c[0]] + [c[i] for i in range(1, l)]
            bisect.insort(ordinals_t_l_1, c[t])
            ordinals_t_l = ordinals_t_l_1
            assert (len(ordinals_t_l)-1) == l

            # Optimal partition of size l on the first t clumps of D
            P_t_l = GetPartitionFromOrdinals(D, ordinals_t_l_1, axis='x')
            I[t][l] = H(Q) + H(P_t_l) - H(P_t_l, Q)

    for l in range(k+1, x+1): I[k][l] = I[k][k]
    
    return I[k][2:x+1]
                  
def GetPartitionOrdinals(D, P, axis='x'):
    P_tilde = GroupPointsByPartition(P)
    if axis == 'x':
        return [-1] + [D.index(get_rightest_point(P_tilde[k])) for k in sorted(P_tilde.keys())]
    elif axis == 'y':
        return [-1] + [D.index(get_uppest_point(P_tilde[k])) for k in sorted(P_tilde.keys())]
  
def GetPartitionFromOrdinals(D, ordinals, axis='x'):
    if not is_sorted_increasing_by(D, 'x'): D = sort_D_increasing_by(D, 'x')
    P = {}
    
    current_partition = 0
    for i, j in pairwise(ordinals):
        from_point = i + 1
        to_point = j
        for p_index in range(from_point, to_point + 1):
            P[D[p_index]] = current_partition
        
        current_partition += 1
    
    return P

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

'''
Computations
'''
def H(P, *Q):

    if not Q:
        return entropy(GetProbabilityDistribution(P)) if isinstance(P, Mapping) else entropy(P)   
    else:
        # Not implemented for len(Q) > 1
        assert len(Q) == 1
        
        Q = Q[0]
        
        # Number of total points
        n_points = len(set(P.iterkeys()) | set(Q.iterkeys()))
        
        # Probability vector for the P-by-Q grid
        G = GetGridMatrix(P, Q).flatten() 
        probabilities = G / float(n_points)
        
        return entropy(probabilities)

def entropy(probs): 
    return -sum(p * log(p, 2) for p in probs if p > 0)
                 
def I(P, Q):
    return H(P) + H(Q) - H(P, Q)

def GetProbabilityDistribution(P):
    partittions = GroupPointsByPartition(P)
    #The probability mass of each partition is the fraction of points that lie in this partition
    prob_mass = lambda p: len(partittions[p]) / float(len(P))
    return map(prob_mass, partittions)

'''
Utils
'''
p_x, p_y = lambda p: p[0], lambda p: p[1]
def get_rightest_point(points): return max(points, key=p_x)
def get_uppest_point(points): return max(points, key=p_y)
def last_abscissa(x_bin): return p_x(get_rightest_point(x_bin))
def last_ordinate(y_bin): return p_y(get_uppest_point(y_bin))

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

def sort_D_increasing_by(D, increasing_by='x'):
    assert increasing_by == 'x' or increasing_by == 'y'
    
    return sorted(D, key=p_x) if increasing_by == 'x' else sorted(D, key=p_y)


def GetGridMatrix(P, Q):
    """
    Each matrix element equals the number of points in the corresponding grid cell.
    """
    P = GroupPointsByPartition(P)
    Q = GroupPointsByPartition(Q)
   
    grid_matrix = np.zeros(shape=(len(Q.keys()), len(P.keys())), dtype=int)
    num_rows, num_columns = grid_matrix.shape[0], grid_matrix.shape[1]
    for r in range(num_rows):
        for c in range(num_columns):
            grid_matrix[r][c] = len(set(Q[r]) & set(P[c]))
    flipped = np.flipud(grid_matrix)
    return flipped

'''
I/O
'''
def visualize(x_axis_parition={}, y_axis_partition={}, step=0.2):
    points = set(chain(x_axis_parition.iterkeys(), y_axis_partition.iterkeys()))
    
    x_axis_parition = GroupPointsByPartition(x_axis_parition)
    y_axis_partition = GroupPointsByPartition(y_axis_partition)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    #Scatter points
    ax.scatter(map(p_x, points), map(p_y, points))
    
    x_bin_edge = lambda x_bin :last_abscissa(x_bin) + step
    y_bin_edge = lambda y_bin:  last_ordinate(y_bin) + step
    
    x_ticks = map(x_bin_edge, x_axis_parition.itervalues())
    y_ticks = map(y_bin_edge, y_axis_partition.itervalues())
    
    ax.get_xaxis().set_ticks(x_ticks)
    ax.get_yaxis().set_ticks(y_ticks)
    
    #Format grid appearance
    ax.grid(True,  alpha=0.5, color='red', linestyle='-', linewidth=1.5)
    
    x_partition_size = len(x_axis_parition.values())
    y_partition_size = len(y_axis_partition.values())
    plt.title(str(x_partition_size)+' - by - '+str(y_partition_size) + ' Grid')
    plt.show()