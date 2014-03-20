from utils import *
from mine import *
from matplotlib import cm

def plot_example_grid():
    x = np.arange(20)
    y = x ** 2 - 1
    m = MINE(x, y)
    q = m.equipartition_y_axis(m.Dy, 3)
    p, _ = m.get_clumps_partition(q)
    plot_partitions(p, q, 'example_grid.png')

def plot_char_matrix_surface(m):
    x = np.arange(m.shape[0])
    y = np.arange(m.shape[1])
    xx, yy = np.meshgrid(x,y)

    @np.vectorize
    def char_value(x,y):
        return m[x][y]


    z = char_value(xx,yy)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    surf = ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
    ax.set_zlim3d(-1.01, 1.01)

    fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.show()


if __name__ == "__main__":
    plot_example_grid()