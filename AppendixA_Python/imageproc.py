import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def get_image_max(file,n):
    A = plt.imread(file)
    A = np.mean(A[:,:,0:2],2)
    A = gaussian_filter(A, np.round(n/2)) 
    maxval = np.max(A) 
    return maxval


