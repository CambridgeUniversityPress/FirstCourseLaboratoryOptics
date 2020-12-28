# axial_irradiance.py
# Author: A. Gretarsson
#
# Plots the maximum irradiance (arb. units) in a beam versus propagation
# distance. Assumes the image filenames have purely numerical names corresponding to
# the location of the image taken (in whatever units used). An attempt is made
# to reduce laser speckle by smoothing the image. This requires the user to
# be judicious in the choice of nsmooth, the linear smoothing size in pixels. 
# 
# Requires: get_image_max

import numpy as np
import matplotlib.pyplot as plt
from imageproc import get_image_max
import os



image_folder = 'sample_images'   # relative or absolute path to the images directory
image_extension = '.tif'         # filename extension ofthe images
nsmooth = 32                     # number of pixels over which to smooth the images

images = os.listdir(image_folder)
posvals = np.array([])
maxvals = np.array([])

for s in range(1,len(images)):
    
    fdir,fnameext = os.path.split(images[s])
    fname,fext = os.path.splitext(fnameext)
    fullpath = os.path.join(image_folder, fnameext)
    if fext==image_extension:
         posvals = np.concatenate((posvals,\
                   [float(fname)]));
         maxvals =  np.concatenate((maxvals,\
                   [get_image_max(fullpath,nsmooth)]))

plt.plot(posvals,maxvals,'bs',linewidth=2);
plt.grid(True)
plt.xlabel('Position  ( mm )');
plt.ylabel('Max. Irradiance  ( arb. units )');
plt.show()