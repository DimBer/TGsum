import numpy as np 


E = np.genfromtxt('graph', dtype=int)

N = E.shape[0]

for i in range(9):
    p = int(((i+1)/10)*N)
    s = np.random.choice( range(N), size=p, replace=False)
    np.savetxt('samples/'+str(i+1), E[s,:], fmt='%d %d')


