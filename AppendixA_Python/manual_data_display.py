# Demonstrates how to enter and display data
import numpy as np
import matplotlib.pyplot as plt

data = np.array(# \=line continnuation char.\
 [ #  X       deltaX  Y       deltaY \
  [   0.9,    0.2,    2.1,    0.2 ],\
  [   1.4,    0.3,    3.4,    0.3 ],\
  [   2.1,    0.2,    3.1,    0.6 ],\
  [   2.7,    0.2,    4.8,    0.5 ],\
  [   3.4,    0.2,    5.1,    0.2 ]\
 ])
    
deltaX = data[:,1]         # X errors
deltaY = data[:,3]         # Y errors  

plt.errorbar(data[:,0],data[:,2],\
    xerr=deltaX,yerr=deltaY,fmt='o')
plt.axis([0,4,0,7])
plt.grid(axis='both')

plt.plot([0.5,3.75],[1.5, 5.8],\
         'r-');            # "by eye" fit
plt.xlabel('$\Delta L$   ( nm )')
plt.ylabel('$T$   ( $^\circ C$ )')
plt.show()