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


from collections import Counter

import numpy as np
from utils import *


def H(distribution):
    def entropy(P):
        '''
        Return the Shannon entropy of a probability vector P
        See http://www.scholarpedia.org/article/Entropy#Shannon_entropy
        '''
        h = -np.fromiter((i * np.log2(i) for i in P if i > 0), dtype=np.float64).sum()
        return h
    return entropy(distribution.ravel())

def HQ(Q):
    n = len(Q)
    histogram = np.fromiter(Counter(Q.itervalues()).itervalues(), dtype=np.float64)
    print histogram
    return H(histogram / np.float64(n))

def HPQ(P, Q):
    return H(GetGridHistogram(P, Q))

def HP(P):
    return H(get_partition_histogram(P))

def getXDistribution(grid_histogram):
    return grid_histogram.sum(axis=0)

def getYDistribution(grid_histogram):
    return grid_histogram.sum(axis=1)

def I(joint_distribution_histogram):
    x_distribution = getXDistribution(joint_distribution_histogram)
    y_distribution = getYDistribution(joint_distribution_histogram)
    joint_distribution = joint_distribution_histogram.ravel()
    return H(x_distribution) + H(y_distribution) - H(joint_distribution)

