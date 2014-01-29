import matplotlib.pyplot as plt
import numpy as np
from utils import *
from mine import *

def plot_example_grid():
    #Run, visualize and export an example grid
    x = np.linspace(-100, 100, 100)
    y = x**2 + x
    D = zip(x, y)
    Dx = sort_D_increasing_by(D, 'x')
    Dy = sort_D_increasing_by(D, 'y')
    Q = EquipartitionYAxis(Dy, y=10)
    
    P = GetClumpsPartition(Dx, Q)
    x_axis_parition = GroupPointsByPartition(P)
    y_axis_partition = GroupPointsByPartition(Q)
    points = set()
    
    for partition in x_axis_parition.values():
        for p in partition:
            points.add(p)
            
    for partition in y_axis_partition.values():
        for p in partition:
            points.add(p)

    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # Scatter points
    ax.scatter(map(p_x, points), map(p_y, points))
    
    x_bin_edge = lambda x_bin :last_abscissa(x_bin) + 0.2
    y_bin_edge = lambda y_bin:  last_ordinate(y_bin) + 0.2
    
    x_ticks = map(x_bin_edge, x_axis_parition.values())
    y_ticks = map(y_bin_edge, y_axis_partition.values())
    
    ax.get_xaxis().set_ticks(x_ticks)
    ax.get_yaxis().set_ticks(y_ticks)
    
    # Format grid appearance
    ax.grid(True, alpha=0.5, color='red', linestyle='-', linewidth=1.5)
    
    x_partition_size = len(x_axis_parition.values())
    y_partition_size = len(y_axis_partition.values())
    plt.title(str(x_partition_size) + ' - by - ' + str(y_partition_size) + ' Grid')
    #plt.show()
    plt.savefig("../doc/example_grid.png")
    

plot_example_grid()