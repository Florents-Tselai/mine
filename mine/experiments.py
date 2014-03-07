from encodings.cp437 import encoding_map
from collections import OrderedDict
from time import time

from numpy.core.fromnumeric import size
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


"+chr(94)+"
functions = OrderedDict()
functions['Linear'] = "$y=x$"
functions['Parabolic'] = "$y=4*(x-1/2)^2$"
functions['Cubic'] = "$128(x-1/3)^2-48(x-1/3)" + chr(94) + "2-12(x-1/3)+2$"
functions['Exponential'] = "$y=10^(10x)-1$"
functions['Linear/Periodic'] = "$y=\sin(10\pi x)+x$"
functions['Sinusodial (Fourier Frequency)'] = "$y=\sin(16\pi x)$"
functions['Sinusodial (non-Fourier Frequency)'] = "$y=sin(13\pi x)$"
functions['Sinusodial (Varying Frequency)'] = "$y=\sin(7\pi x(1+x)$"
functions['Categorical'] = "$64 points chosen from: {}$"
functions['Random'] = "$random number generator$"


def get_y(x, f):
    if f == 'Linear':
        return x;
    if f == 'Parabolic':
        return 4 * (x - 1. / 2) ** 2
    if f == 'Cubic':
        return 128 * (x - 1. / 3) ** 3 - 48 * (x - 1. / 3) - 12 * (x - 1. / 3) + 2
    if f == 'Exponential':
        return 10 ** (10 * x) - 1
    if f == 'Linear/Periodic':
        return np.sin(10 * np.pi * x) + x
    if f == 'Sinusodial (Fourier Frequency)':
        return np.sin(16 * np.pi * x)
    if f == 'Sinusodial (non-Fourier Frequency)':
        return np.sin(13 * np.pi * x)
    if f == 'Sinusodial (Varying Frequency)':
        return np.sin(7 * np.pi * x * (1 + x))
    if f == 'Categorical':
        return x
    if f == 'Random':
        return x


def run(f, n):
    x = np.linspace(0, 1, num=n)
    y = get_y(x, f)
    t0 = time()
    from mine import MINE

    m = MINE(x, y)
    b = pow(len(x), 0.6)
    M = m.approx_char_matrix(m.D, b)
    mic = np.max(M[~np.isnan(M)])

    t = time()
    return np.round(t - t0, 3)


def experiment_1():
    n_points = [100, 200, 300, 400, 500]
    func_df = pd.DataFrame(index=functions.keys())
    df = pd.DataFrame(index=functions.iterkeys(), columns=n_points)

    for function, exp in df.iterrows():
        for n in df.columns:
            exp[n] = run(function, n / 5)

    df.T.plot(sharex=True)

    plt.xlabel("Number of points")
    plt.ylabel("Running time")
    plt.legend(loc='upper left', prop={'size': 10})
    plt.savefig("/home/florents/workspace/mine/doc/experiments/experiment_1_plot.png", dpi=1000)

    df.insert(0, 'Relationship Name', df.index)

    output = df.to_latex(index=False)
    output = output.replace("\$", "$").replace("\\textasciicircumx", "").replace("\\textasciicircum", "")
    text_file = open("/home/florents/workspace/mine/doc/experiments/experiment_1.tex", "w")
    text_file.write(output)
    text_file.close()


def print_functions_definition():
    df = pd.DataFrame();
    df['Relationship Name'] = functions.keys()
    df['Description'] = df['Relationship Name'].map(functions.get)
    output = df.to_latex(index=False)
    output = output.replace("\$", "$").replace("\\textasciicircumx", "").replace("\\textasciicircum", "").replace(
        "\\textbackslash", '\\')
    text_file = open("/home/florents/workspace/mine/doc/experiments/functions_definitions.tex", "w")
    text_file.write(output)
    text_file.close()


if __name__ == "__main__":
    print_functions_definition()
    experiment_1()