from collections import defaultdict, Mapping
import matplotlib.pyplot as plt
from itertools import chain
import numpy as np
from math import log
from copy import copy
from pprint import pprint

"""
Each point is a tuple, thus:
P = (x, y) ==> 
P[0] = x and P[1] = y
"""

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
    if x_ticks: del[x_ticks[len(x_ticks) - 1]]
    if y_ticks: del[y_ticks[len(y_ticks) - 1]]
    
    
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
        S = [p for p in D if p[1] == D[i][1]]
        
        lhs = abs(float(sharp) + float(len(S)) - desiredRowSize)
        rhs = abs(float(sharp) - desiredRowSize)
        
        if ((sharp != 0) and (lhs >= rhs)):
            currRow = currRow + 1
            sharp = 0
            temp1 = float(n) - float(i)
            temp2 = float(y) - float(currRow)
            desiredRowSize = temp1 / temp2
        
        for j in range(0, len(S)): Q[D[i + j]] = currRow + 1
        
        i += len(S)
        sharp += len(S)
    
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
            for j in range(0, s):
                Q_tilde[D[i + j]] = c
            c -= 1
        i = i + s
    
    i = 0
    P = {}
    P[D[0]] = 0 + 1
    for j in range(1, n):
        if Q_tilde[D[j]] != Q_tilde[D[j - 1]]:
            i = i + 1
        P[D[j]] = i + 1
    
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

def H(P=None, Q=None):
    assert P is not None or Q is not None
    
    # TODO: Refactor code duplication
    if P is not None and Q is None:
        if isinstance(P, Mapping):
            return entropy(GetProbabilityDistribution(P))
        else:
            return entropy(P)
            
    elif Q is not None and P is None:
        if isinstance(Q, Mapping):
            return entropy(GetProbabilityDistribution(P))
        else:
            return entropy(Q)
    
    else:
        assert P is not None and Q is not None
        # Compute joint entropy
        n_points = float(len(set(chain(P.iterkeys(), Q.iterkeys()))))
        
        grid_matrix = GetGridMatrix(P, Q)
        probabilities = grid_matrix.flatten() / n_points
        
        return entropy(probabilities)

def entropy(probs): 
    return -sum(p * log(p, 2) for p in probs if p > 0)
                 
def I(P, Q):
    return H(P=P) + H(Q=Q) + H(P=P, Q=Q)

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
    return d

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
    for r in range(0, grid_matrix.shape[0]):
        for c in range(0, grid_matrix.shape[1]):
            grid_matrix[r][c] = len(set.intersection(set(Q[r + 1]), set(P[c + 1])))
    flipped = np.flipud(grid_matrix)
    return flipped
            
