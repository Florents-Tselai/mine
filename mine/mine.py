from collections import defaultdict
import matplotlib.pyplot as plt

import numpy as np
from math import log

"""
Each point is a tuple, thus:
P = (x, y) ==> 
P[0] = x and P[1] = y
"""

p_x, p_y = lambda p: p[0], lambda p: p[1]

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
    
    #Draw partition lines for each axis
    ax.get_xaxis().set_ticks([D_sorted_by_x[c][0] + step for c in x_partition_endpoint_indices])
    ax.get_yaxis().set_ticks([D_sorted_by_y[c][1] + step for c in y_partition_endpoint_indices])
   
    ax.grid(True)
    
    plt.show()
    
def visualize_partition_by_endpoints(D, x_endpoints, y_endpoints, step=0.2):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter([p[0] for p in D], [p[1] for p in D])
    
    #Draw partition lines for each axis
    ax.get_xaxis().set_ticks([p[0]+step for p in x_endpoints])
    ax.get_yaxis().set_ticks([p[1]+step for p in y_endpoints])
   
    ax.grid(True)
    
    plt.show()
    
    pass

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
    
    i = 0
    c = -1 
    
    while(i < n):
        s = 1
        flag = False
        for j in range(i + 1, n):
            if D[i][0] == D[j][0]:
                s += 1
                if Q[D[i]] != Q[D[j]]:
                    flag = True
            else:
                break
            
        if s > 1 and flag:
            for j in range(0, s):
                Q[D[i + j]] = c
            c -= 1
        i = i + s
    
    i = 0
    P = {}
    P[D[0]] = 0 + 1
    for j in range(1, n):
        if Q[D[j]] != Q[D[j - 1]]:
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
    
    if P is not None:
        assert Q is None
        return -sum(p * log(p, 2) for p in P)
    
    elif Q is not None:
        assert P is None
        return -sum(q * log(q, 2) for q in Q)
    
    else:
        assert P is not None and Q is not None
        # Compute joint entropy
        P = np.array(P)
        Q = np.array(Q)
        h = 0.0
        for i in set(P):
            for j in set(Q):
                ppq = np.mean(np.logical_and(P == i, Q == j))
                if ppq > 0:
                    h += ppq * np.log2(ppq)
                else:
                    h += 0.0
        return -h
                    
def I(P, Q):
    P = np.array(P)
    Q = np.array(Q)
    Hpq = H(P, Q)
    h = 0.0
    for j in set(Q):
        for i in set(P):
            ppq = np.mean(np.logical_and(P == i, Q == j))
            pp = np.mean(P == i)
            pq = np.mean(Q == j)
            if ppq > 0 and pp > 0 and pq > 0:
                h += ppq * np.log2(pp * pq)
            else:
                h += 0
    return (-h - Hpq)

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

def GetProbabilityDistribution(P, n):
    """
    n: number of total points
    """
    d = GroupPartitionsPoints(P)    
    return [float(len(d[k])) / float(n) for k in sorted(d.keys())]
