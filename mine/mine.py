from copy import copy
from collections import OrderedDict, defaultdict
import numpy as np
from math import log

def visualize_partition(points, x_partition=[], y_partition=[]):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter([p[0] for p in points], [p[1] for p in points])
    ax.get_xaxis().set_ticks(x_partition)
    ax.get_yaxis().set_ticks(y_partition)
   
    ax.grid(True)
    
    plt.show()
    
def EquipartitionYAxis(D, y):
    n= len(D)
    
    desiredRowSize = float(n) / float(y)
    
    i = 0
    sharp = 0
    currRow = 0
    
    Q = {}
    while(i < n):
        S = [p for p in D if p[1] == D[i][1]]
        
        temp1 = abs(float(sharp) + float(len(S)) - desiredRowSize)
        temp2 = abs(float(sharp) - desiredRowSize)
        
        if ((sharp != 0) and (temp1 >= temp2)):
            currRow = currRow + 1
            sharp = 0
            temp1 = float(n) - float(i)
            temp2 = float(y) - float(currRow)
            desiredRowSize = temp1 / temp2
        
        for j in range(0, len(S)): Q[D[i+j]] = currRow + 1
        
        i += len(S)
        sharp += len(S)
    
    return Q

def GetClumpsPartition(D, Q):
    n = len(D)
    
    i = 0
    c = -1 
    
    while(i < n):
        s = 1
        flag = False
        for j in range(i+1, n):
            if D[i][0] == D[j][0]:
                s += 1
                if Q[D[i]] != Q[D[j]]:
                    flag = True
            else:
                break
            
        if s > 1 and flag:
            for j in range(0, s):
                Q[D[i+j]] = c
            c -= 1
        i = i + s
    
    i = 0
    P = {}
    P[D[0]] = 0 + 1
    for j in range(1, n):
        if Q[D[j]] != Q[D[j-1]]:
            i = i + 1
        P[D[j]] = i + 1
    
    return P

def GetSuperclumpsPartition(D, Q, k_hat):
    pass

def H(P=None, Q=None):
    if P is not None:
        assert Q is None
        return - sum(p*log(p,2) for p in P)
    elif Q is not None:
        assert P is None
        return - sum(q*log(q,2) for q in Q)
    else:
        assert P is not None and Q is not None
        #Compute joint entropy
        probs = []
        for p in set(P):
            for q in set(Q):
                probs.append(np.mean(np.logical_and(P == p, Q == q)))
                return np.sum(-p * np.log2(p) for p in probs)

def getPartitionIndices(P, axis='x', step=0.3):
    d = GetPartitionGroups(P)
    
    indices = []
    for k in sorted(d.keys()):
        if axis == 'x':
            indices.append(max([p[0] for p in d[k]]) + step)
        elif axis == 'y':
            indices.append(max([p[1] for p in d[k]]) + step)
            
    del indices[len(indices)-1]
    
    return indices

def GetPartitionGroups(P):
    d = defaultdict(list)
    for k,v in P.iteritems(): 
        d[v].append(k)
    return d

def GetProbabilityDistribution(P, n):
    """
    n: number of total points
    """
    d = GetPartitionGroups(P)    
    return [float( len(d[k])) / float(n) for k in sorted(d.keys())]

D = [(1,1), (1,2), (1,3), (1,4), (2,3), (2,4), (3,5), (4,6), (5,6), (6,6), (7,5), (8,3), (9,2), (9,1)]
x_partition = [1.8, 2.2, 7.8, 8.2]
y_partition = [2.5, 4.6]