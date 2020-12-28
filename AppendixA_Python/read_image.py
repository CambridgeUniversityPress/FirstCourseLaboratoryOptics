import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

A = mpimg.imread('myphoto.jpg');

A=np.mean(A,2)
plt.pcolormesh(A,\
	 shading='flat',cmap='bone');