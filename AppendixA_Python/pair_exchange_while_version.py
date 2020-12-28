import numpy as np

vec = np.arange(1,10)
s = 1
while s <= len(vec)-1:
    tmp = vec[s]
    vec[s] = vec[s-1]
    vec[s-1] = tmp
    s=s+2
    
print(vec)


