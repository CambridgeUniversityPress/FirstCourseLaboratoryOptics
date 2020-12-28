import numpy as np

vec = np.arange(1,10);
for s in range(1,len(vec),2):
    tmp = vec[s]
    vec[s] = vec[s-1]
    vec[s-1] = tmp
    
print(vec)


