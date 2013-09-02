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
    D = sorted(D, key=lambda p: p[1])
    print D
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
                if Q[D[i]] == Q[D[j]]:
                    flag = True
                    s += 1
            else:
                break
            
        if s > 1 and flag == True:
            for j in range(0, s):
                Q[D[i+j]] = c
                c -= 1
        i += s
    
    i = 0
    P = {}
    P[D[0]] = 0
    for j in range(1, n):
        if Q[D[j]] != Q[D[j-1]]:
            i += 1
        P[D[j]] = i
        
    return P
                    

    
D = [(1,1), (1,2), (1,3), (1,4), (2,3), (2,4), (3,5), (4,6), (5,6), (6,6), (7,5), (8,3), (9,2), (9,1)]
x_partition = [1.8, 2.2, 7.8, 8.2]
y_partition = [2.5, 4.6]