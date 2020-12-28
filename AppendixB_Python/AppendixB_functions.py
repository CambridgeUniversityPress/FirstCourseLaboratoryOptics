# AppendixB_functions.py
# Author: A. Gretarsson
#
# This file contains the Python functions for Appendix B in
# "A First Course in Laboratory Optics" by A. Gretarsson. These are 
# Python versions of the Matlab functions shown there.

import numpy as np
import matplotlib.pyplot as plt

def beamradius(params,z):
    """Returns the field radius of a TEM_00 mode beam at any point z 
    along the optic axis. 
    
    SYNTAX: w,R,zR = beamradius([w0,zw,lam],z);
    
    w0 = waist size
    zw = position of waist
    lam = wavelength
    
    w  = spot size (field radius) at z
    R  = curvature of phasefront at z
    zR = Raleigh length.
    
    The input arguments w0, zw, lam, z all need to be in the same 
    units. The output arguments will be in those units."""
    
    w0=params[0]                        # beam width at waist [e.g. meters]
    zw=params[1]                        # waist position [e.g. meters]
    lam = params[2]                     # wavelength [meters]
    
    zR=np.pi*w0**2/lam                  # Raleigh length [e.g. meters]
    w=w0*np.sqrt(1+((z-zw)/zR)**2)      # beam width at z [e.g. meters]
    R=z*(1+(zR/z)**2)                   # beam phasefront curvature at z

    return  w,R,zR                      # values at pos z [e.g. meters]

def prop(q1,abcd,mode=[0,0],p1=1):
    """Propagates a Gaussian beam (TEM_nm) with complex radius of 
    curvature q1 and amplitude factor p1 (optional), using the abcd 
    matrix supplied.
    
    Returns the new complex beam radius q=(A*q1+B)/(C*q1+D) and the 
    new amplitude factor p = 1/(A+B/q1)^(1+n+m) by which the field is
    multiplied. If q1 is a vector q and p will be vectors of the same size.

    SYNTAX: [q,p]=prop(q1,abcd <,[n,m],p1>);
            <...> indicates optional arguments
            
    For a Hermite Gaussian n,m are the mode designators."""

    A=abcd[0][0]
    B=abcd[0][1]
    C=abcd[1][0]
    D=abcd[1][1]
    
    n=mode[0]
    m=mode[1]
    
    q = (A*q1 + B)/(C*q1 + D)
    p = p1*np.exp(1j*np.angle(1/(A+B/q1)**(1+n+m)))
    
    return q,p

def q_(w,R,lam=1064.0e-9):
    """Returns the q-factor of a Gaussian beam given the spot size, w,
    phasefront radius of curvature, R, and wavelength, lam.

    SYNTAX: q=q_(w,R <,lam>);
                    <...> indicates optional arguments

    w     = 1/e Field radius 
    R     = Radius of curvature of phasefront
    lam   = wavelength 

    w, R and lam must all be in the same units."""

    if R!=np.inf:
        q=np.pi*w**2*R/(np.pi*w**2-1j*R*lam)
    else:
        q=1j*np.pi*w**2/lam

    return q

def R_(q,lam=1064.0e-9):
    """﻿Returns the phasefront radius of curvature, R, and the beam 
    width, w, of a Gaussian beam. Accepts the complex beam radius,
    q, and the wavelength. 

    SYNTAX: R,w=R_(q <,lambda>);   
                <...> indicates optional arguments

    q       = q-factor of the beam at the position where R and w are to
              be found. q can be a vector
    lam     = wavelength. Can be a vector or scalar.
    w       = beam radius
    R       = beam phasefront curvature"""
    
    if np.real(q)!=0: # provided q isn't purely imaginary, use usual formulas
        w = np.sqrt(lam/np.pi * np.imag(q)*(1+np.real(q)**2/np.imag(q)**2))
        R=np.real(q)*(1+np.imag(q)**2/np.real(q)**2)*np.ones(np.shape(w))
    else:             # if q is purely imaginary, R will be infinite
        w = np.sqrt(lam/np.pi * np.imag(q)*(1+np.real(q)**2/np.imag(q)**2))
        R = np.inf
        
    return R,w
    
def imdespeckle(imagefile, threshold):
    """Performs a 2D fourier transform on the image in "imagefile", then sets to zero, all spatial 
    frequency bins with amplitude below the value given in  "threshold". The result is then 
    fourier transformed back to the spatial domain to produce the despeckled image. The beam 
    occupies a small part of the freq. domain image; speckle occupies the rest. Speckle has low 
    power per frequency bin. Suppressing all bins with low power therefore gets rid of speckle.

    SYNTAX:  I = imdespeckle(imagefile, threshold);

    NOTES: - threshold = 1 is usually a good starting point.
           - The image is assumed to be a monochrome image. If it's not, it's converted
             to a monochrome image before processing. 
           
    INPUT ARGUMENTS
    ---------------
    imagefile:    Full path and filename to the image (any format readable by "imread").
    threshold:    Frequency components with log10(mag) below threshold are discarded."""


    data = plt.imread(imagefile);           # image is read into the array "data"
    data = np.mean(data,2);                 # convert to greyscale
    
    # Perform the 2D numerical fourier transform and scale it correctly. The result is a
    # picture of the image in "frequency space" (spatial frequency, that is).
    N1 = np.shape(data)[0]                              # number of rows
    N2 = np.shape(data)[1]                              # number of columns
    F=np.fft.fftshift(np.fft.fft2(data)/np.sqrt(N1*N2)) # 2D FT with zero freq's in center

    # Threshold the fourier transformed image
    pixels_below_threshold = np.log10(np.abs(F))<threshold # logical mask for pixels -> 0
    Fthresh = F                                         # start unthresholded
    Fthresh[pixels_below_threshold] = 0                 # set pixels below threshold to 0                   
    
    # Finally, perform the inverse transform on the thresholded data to get back
    # to position space. (I.e. to get back our image.).
    despekld_image = np.abs(np.fft.ifft2(Fthresh)*np.sqrt(N1*N2))

    # Now display the results
    plt.figure(1)                                       # open figure 1
    plt.clf()                                           # clear it in case previously used
    ax1 = plt.axes()                                    # define a set of axes
    ax1.pcolormesh(despekld_image, cmap='bone')         # plot the despeckled image
    ax1.set_aspect('equal', 'box')                      # set aspect ratio to be correct
    ax1.set_title('Despeckled Image')                   # add a title
    plt.show()                                          # display the plot

    fig2 = plt.figure(2)
    plt.clf()
    ax2 = plt.axes()
    with np.errstate(divide='ignore'):                  # suppresses warning for "log10(0)" 
        c2 = ax2.pcolormesh(np.log10(np.abs(Fthresh)), cmap='viridis') # plot the FT
    fig2.colorbar(c2)
    ax2.set_aspect('equal', 'box')
    ax2.set_title('Log10 of the 2D FFT, Thresholded')
    plt.show()
   
    return despekld_image
    