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
    D = sorted(D, key=lambda p: p[1], reverse=True)
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
        
        for j in range(0, len(S)): Q[D[i+j]] = currRow
        
        i += len(S)
        sharp += len(S)
    
    return Q

def GetClumpsPartition(D, Q):
    pass


    
D = [(1,1), (1,2), (1,3), (1,4), (2,3), (2,4), (3,5), (4,6), (5,6), (6,6), (7,5), (8,3), (9,2), (9,1)]
x_partition = [1.8, 2.2, 7.8, 8.2]
y_partition = [2.5, 4.6]

for k,v in EquipartitionYAxis(D, 3).iteritems():
    print k,v