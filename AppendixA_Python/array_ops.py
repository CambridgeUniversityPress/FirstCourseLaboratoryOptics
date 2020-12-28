# This code is intended to be added line by line to the console/command window
import numpy as np

A = np.array([ [2,0,1], [0,1,3],  [-1,1,0] ])
B = np.array([ [1,1,2], [1,2,-2], [1,4,1] ])

A[1,0]      # Elem. at second row first col.
A[1,:]      # Entire second row
A[:,1]      # Entire second column
A.T         # Transpose
A @ B       # Matrix multiplication
A * B       # Element-wise mult.