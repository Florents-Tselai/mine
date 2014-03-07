from utils import *
from mine import *


def plot_example_grid():
    x = np.arange(20)
    y = x ** 2 - 1
    m = MINE(x, y)
    q = m.equipartition_y_axis(m.Dy, 3)
    p, _ = m.get_clumps_partition(q)
    plot_partitions(p, q, 'example_grid.png')


if __name__ == "__main__":
    plot_example_grid()