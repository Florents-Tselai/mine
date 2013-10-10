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
from mine_utils import *
from mine_computations import *
import matplotlib.pyplot as plt
import numpy as np


def ApproxMaxMI(D, x, y, k_hat):
    assert x > 1 and y > 1 and k_hat > 1 
    
    Q = EquipartitionYAxis(D, y)
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
    c = GetPartitionOrdinals(D, super_clumps_partition, axis='x')
    
    # Total number of clumps
    k = len(set(super_clumps_partition.values()))
    assert k == len(c) - 1

    # Find the optimal partitions of size 2 
    I = np.zeros(shape=(k + 1, x + 1))
    
    for t in xrange(2, k + 1):
        s = max(xrange(1, t + 1),
                key=lambda s: 
                            H(GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]], axis='x')) - 
                            H(GetPartitionFromOrdinals(D, ordinals=[c[0], c[s], c[t]], axis='x'), Q))
        
        # Optimal partition of size 2 on the first t clumps
        optimal_2_partition_ordinals = [c[0], c[s], c[t]]
        P_t_2 = GetPartitionFromOrdinals(D, ordinals=optimal_2_partition_ordinals, axis='x')
        I[t][2] = H(Q) + H(P_t_2) - H(P_t_2, Q)
   
    # Inductively build the rest of the table of optimal partitions
    for l in xrange(3, x + 1):
        for t in xrange(l, k + 1):
            
            s = max(range(l - 1, t + 1),
                    key=lambda s: float((c[s] / c[t])) * (I[s][l - 1] - H(Q)) - float(((c[t] - c[s]) / c[t])) * H(GetPartitionFromOrdinals(D, [c[s], c[t]], axis='x'), Q)
                    )
            ordinals_t_l_1 = [c[0]] + [c[i] for i in xrange(1, l)]
            bisect.insort(ordinals_t_l_1, c[t])
            ordinals_t_l = ordinals_t_l_1
            assert (len(ordinals_t_l) - 1) == l

            # Optimal partition of size l on the first t clumps of D
            P_t_l = GetPartitionFromOrdinals(D, ordinals_t_l_1, axis='x')
            I[t][l] = H(Q) + H(P_t_l) - H(P_t_l, Q)

    for l in xrange(k + 1, x + 1): I[k][l] = I[k][k]
    
    return I[k][2:x + 1]
