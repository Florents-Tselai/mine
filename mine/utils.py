p_x, p_y = lambda p: p[0], lambda p: p[1]

def get_rightest_point(points): 
    return max(points, key=p_x)

def get_leftest_point(points): 
    return min(points, key=p_x)

def get_uppest_point(points): 
    return max(points, key=p_y)

def get_downest_point(points): 
    return min(points, key=p_y)

def last_abscissa(x_bin): 
    return p_x(get_rightest_point(x_bin))

def last_ordinate(y_bin): 
    return p_y(get_uppest_point(y_bin))

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