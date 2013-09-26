from collections import defaultdict, Mapping
import matplotlib.pyplot as plt
from itertools import chain
import numpy as np
from math import log
from copy import copy
from pprint import pprint
import bisect

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
        s = len([p for p in D if p[1] == D[i][1]])
        
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

    c = GetPartitionOrdinals(GetSuperclumpsPartition(D, Q, k_hat), D, axis='x')
    
    #Find the optimal partitions of size 2
    k = len(c) - 1
    I = np.array(shape=(k, x))
    for t in range(2, k+1):
        s = max(range(1,t+1), 
                key=lambda s: H(P=GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]])) - H(P=GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]]), Q=Q)
                )
        P_t_2 = GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]])
        I[t][2] = H(Q=Q) + H(P=P_t_2) - H(P=P_t_2, Q=Q)
    
    #Inductively build the rest of the table of optimal partitions
    for l in range(3, x+1):
        for t in range(l, k+1):
            s = max(range(l-1, t+1), 
                    key=lambda s: float((c[s]/c[t])) * (I[s][l-1] - H(Q)) - float(((c[t]-c[s]) / c[t])) * H([c[s], c[t]], Q)
                    )
            ordinals_s_l_1 = [s, l-1]
            bisect.insort(ordinals_s_l_1, c[t])
            P_t_l = GetPartitionFromOrdinals(D, ordinals_s_l_1)
            I[t][l] = H(Q=Q) + H(P=P_t_l) - H(P=P_t_l, Q=Q)
    
    #for l in range(k+1, x+1): 
    for l in range(k+1, x+1):
        I[k][l] = I[k,k]
    return I
            
        
def GetPartitionOrdinals(D, P, axis='x'):
    #TODO: Fix bug! We need a copy here
    P = GroupPartitionsPoints(P)
    ordinals = [D.index(get_rightest_point(P[k])) for k in sorted(P.keys())]
    #We don't need the last one
    del ordinals[-1]
    return ordinals
    

def GetPartitionFromOrdinals(D, ordinals, axis='x'):
    P = {}
    
    if len(ordinals) == 3:
        for i in range(ordinals[0]+1):
            P[D[i]] = 0
        for i in range(ordinals[0]+1, ordinals[1]+1):
            P[D[i]] = 1
        for i in range(ordinals[1]+1, ordinals[2]+1):
            P[D[i]] = 2
        for i in range(ordinals[2]+1, len(D)):
            P[D[i]] = 3
    if len(ordinals) == 2:
        for i in range(0, ordinals[0]+1):
            P[D[i]] = 0
        for i in range(ordinals[0]+1, ordinals[1]+1):
            P[D[i]] = 1
        for i in range(ordinals[1]+1, len(D)):
            P[D[i]] = 2
    return P
        
    

def Hp3(D, c_0, c_s, c_t):
    
    pass

def Hp3Q(c_0, c_s, c_t, Q):
    pass


def H(P, *Q):

    if not Q:
        return entropy(GetProbabilityDistribution(P)) if isinstance(P, Mapping) else entropy(P)
    
    else:
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
    num_rows, num_columns = grid_matrix.shape[0], grid_matrix.shape[1]
    for r in range(num_rows):
        for c in range(num_columns):
            grid_matrix[r][c] = len(set(Q[r]) & set(P[c]))
    flipped = np.flipud(grid_matrix)
    return flipped
            
