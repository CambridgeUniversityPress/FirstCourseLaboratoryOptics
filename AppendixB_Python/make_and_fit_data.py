# -------------------------------------------------------------------------------------------
# DATA FITTING DEMONSTRATION 
# Language: Matlab (see Line 158 for modifications needed to run under Octave)
# -------------------------------------------------------------------------------------------
# This script makes some simulated data for the intensity at the output of a michelson as a 
# function of mirror position with added noise.  It then fits the data to a straight line
# using one of Matlab's built-in curve-fitting tools: lsqnonlin. It handles non-uniform
# uncertainties correctly and propagates the data uncertainties into the best fit parameters.
# -------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#-------------------------------------------------------------------
# Functions used in script 
#-------------------------------------------------------------------

# This function is to be fit to the data
def fitfunc(x,a1,a2):
    y=a1+a2*x;                              # straight line (can be changed to anything)
    return y


# Convenience function for calculating the chi-square
def ChiSqr(a,fitfunc,x,y,yerr):
    wr=(fitfunc(x,a[0],a[1])-y)/yerr        # chisqr is just the quadrature sum of ...
    C = sum(wr**2);                         # the weighted residuals 
    return C

#-----------------------------------------------------------------
# Input the data
#------------------------------------------------------------------

x = np.linspace(-5,5,9)                     # x is the change in the length of one arm
x = x + 0.2*np.random.randn(np.shape(x)[0]) # Make the vector unevenly spaced

a1_actual = 1/2*(1+(4*np.pi/633)*2)         # Offset (2 nm phase offset)
a2_actual = 2*np.pi/633                     # Slope (nm^-1). Michelson made with a HeNe
y = a1_actual + a2_actual*x                 # The physical relationship betw. x & y

sigma_y = 0.025                             # The std of y values (the actual uncertainty)
y = y + sigma_y*np.random.randn(np.shape(y)[0]) # Add the fluctuations due to the uncertainty
yerr =  np.abs(y)*0.02                      # 2% uncertainty estimate made by the observer

#-----------------------------------------------------------------
# Display the data
#-----------------------------------------------------------------

plt.figure(1);                              # Open a figure in which to plot
plt.clf()                                   # Remove any previous plots
plt.errorbar(x,y,yerr,fmt='none',\
             ecolor=[0.2,0.4,0.8],capsize=3); # Plot the data with the unc. estimates
plt.plot(x,y,'o',mfc='none',color=[0.2,0.4,0.8])
plt.xlabel('Mirror motion (nm)');           # Label the x axis
plt.ylabel('$P_\mathrm{out} / P_\mathrm{one~arm}$'); # Label the y axis
plt.xlim([-10,10]);                         # Choose the x axis limits of the plot
plt.ylim([0.4,0.65]);                       # Choose the y axis limits of the plot]
plt.grid(True)                              # Show grid lines
plt.show()


#-----------------------------------------------------------------
# Fit the data to "fitfunc"                 # see end of script for functions
#-----------------------------------------------------------------
a0 = np.array([0,1]);                       # Initial guess at the best fit values.
a,acov = curve_fit(fitfunc,x,y,p0=a0,sigma=yerr,absolute_sigma=True)   # fit the data
xfit = np.linspace(-8,8,100);
plt.plot(xfit,fitfunc(xfit,a[0],a[1]),'r--')
plt.title('Data and Fit')
residual = y-fitfunc(x,a[0],a[1]);          # Residual is data minus the fit
X2red = 1/(np.shape(x)[0]-np.shape(a)[0])\
    *np.sum((residual)**2./yerr**2);        # Reduced chi-squared ~ 1 if fit is good.
print('')
print('X2red = '+"%0.2f"%X2red);            # Show reduced chisquare in the console

#-----------------------------------------------------------------
# Direct estimate of the parameter uncertainties
#-----------------------------------------------------------------
aerr = np.sqrt(np.diag(acov))      # uncertainties in the parameters
print('')
print('Solution:')
for s in range(len(a)):
    print("%0.4f"%a[s]+' +/- '+"%0.4f"%aerr[s]); # Show best fit values and unc. in console

#-----------------------------------------------------------------
# OPTIONAL:
#
# Make and display chi-square cuts for each fit variable. Also,
# draw the line corresponding to chi-square increasing by 1.
#-----------------------------------------------------------------

a_1 = np.linspace(0.4,0.6,200);             # The domain of a_1-axis chi-square cut
a_2 = np.linspace(0.005,0.015,200)          # The domain of a_2-axis chi-square cut
X2_a1cut = np.zeros(len(a_1));              # will hold the y-values of the a_1 cut
X2_a2cut = np.zeros(len(a_2));              # will hold the y-values of the a_2 cut
for s in range(len(a_1)):
    X2_a1cut[s] = ChiSqr([a_1[s],a[1]],fitfunc,x,y,yerr); # the a_1 cut chi-square values 
for s in range(len(a_1)):
    X2_a2cut[s] = ChiSqr([a[0],a_2[s]],fitfunc,x,y,yerr); # the a_2 cut chi-square values

plt.figure(2);                              # Open a figure in which to plot
plt.clf()                                   # Remove any previous plots
plt.subplot(2,1,1)
plt.plot(a_1,X2_a1cut,'-',a[0],np.min(X2_a1cut),'o',\
    np.array([np.min(a_1),np.max(a_1)]),\
    np.array([np.min(X2_a1cut),np.min(X2_a1cut)])+1,'--') # plot X^2 cut in the a_1 dir, the min, and min+1 line
plt.xlabel('$a_1$');           # Label the x axis
plt.ylabel('$\chi^2$'); # Label the y axis
plt.title('$\chi^2$ Cuts')
plt.grid(True)                              # Show grid lines
plt.show()

plt.subplot(2,1,2)
plt.plot(a_2,X2_a2cut,'-',a[1],np.min(X2_a2cut),'o',\
    np.array([np.min(a_2),np.max(a_2)]),\
    np.array([np.min(X2_a2cut),np.min(X2_a2cut)])+1,'--') # plot X^2 cut in the a_1 dir, the min, and min+1 line
plt.xlabel('$a_2$');           # Label the x axis
plt.ylabel('$\chi^2$'); # Label the y axis
plt.grid(True)                              # Show grid lines
plt.show()

#-----------------------------------------------------------------
# OPTIONAL:
#
# Generate & Display the chi-squared surface
# (This section for illustrative purposes and can be omitted.)
#-----------------------------------------------------------------

# Set up the domain. a1 and a2 are matrixes of parameter values at which to calculate chi-squared
a1,a2 = np.meshgrid(np.linspace(0.47,0.55,100),np.linspace(-0.01,0.03,100));

X2 = np.zeros(np.shape(a1));     # This is the surface we will be finding
for k in range(np.shape(a1)[0]):                     # step through all the values of a1 and a2
    for s in range(np.shape(a2)[1]):                 # in the desired range
        X2[k,s] = ChiSqr([a1[k,s],a2[k,s]],fitfunc,x,y,yerr); # Uses ChiSqr function def. above


fig3 = plt.figure(3)                                # open figure 2
plt.clf()                                           # clear it in case previously used
ax3 = plt.axes()                                    # define a set of axes
X2mindisp = np.min(X2)
X2maxdisp = X2mindisp + 100
c3 = ax3.pcolormesh(a1,a2,X2, cmap='bone',shading='gouraud',\
                    vmin=X2mindisp,vmax=X2maxdisp)  # plot the field plane irradiance
c3b=plt.contour(a1,a2,X2,levels=[np.min(X2)+1],\
    colors=['#FFFFFF'],linestyles='dashed');        # plot one contour at min(X2) + 1
plt.plot(a[0],a[1],'w.')                            # plot a point at the best fit solution
ax3.set_aspect('equal', 'box')                      # set aspect ratio to be correct
ax3.set_title('Field Plane Irradiance')             # add a title
plt.xlabel('$a_1$')                                 # label x axis
plt.ylabel('$a_2$')                                 # label y axis
c3bar = fig3.colorbar(c3)                           # add a color bar
c3bar.set_label('$\chi^2$', rotation=90)            # label the colorbar
plt.show()                                          # display the plot






