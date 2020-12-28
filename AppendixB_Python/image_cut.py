# -------------------------------------------------------------------------------------------
# IMAGE CUT DEMONSTRATION 
# Language: Python 3
# -------------------------------------------------------------------------------------------
# This script loads an image from a file whose full or relative path is specified. The cut
# is made vertically through the center of the image. Then the vector containing
# the cut intensities is smoothed to reduce laser speckle and displayed.
# -------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

imagefile = 'myphoto.jpg';          # User specified path to image file
A=plt.imread(imagefile);            # image is read in as an array A
A = A/255;                          # normalize.
A = np.mean(A[:,:,0:2],2);          # color array has three pages (R,G,B) average->greyscale

col = round(np.shape(A)[1]/2);      # the column index corresponding to the image center
cut = A[:,col];                     # this is how the cut is actually taken

softcut = np.convolve(cut,np.ones(10),'valid')/10 # smooth laser speckle using a moving mean

plt.plot(cut);                      # plot the original data
plt.plot(softcut);                  # and the smoothed version
plt.show()