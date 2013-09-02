from mine import visualize_partition

"""
Example from Albanese et. al. pg. 4
http://mpba.fbk.eu/sites/mpba.fbk.eu/files/albanese12cmine_suppmat.pdf

"""

D = [(1,1), (1,2), (1,3), (1,4), (2,3), (2,4), (3,5), (4,6), (5,6), (6,6), (7,5), (8,3), (9,2), (9,1)]
x_partition = [1.8, 2.2, 7.8, 8.2]
y_partition = [2.5, 4.6]

visualize_partition(D, x_partition, y_partition)

