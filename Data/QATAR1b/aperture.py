'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Test Aperture Determination for Qatar-1 b
'''

from numpy import zeros
from numpy import genfromtxt
from matplotlib.pyplot import *

# Read in data
aperture = genfromtxt('aperture.txt', dtype = 'float')

# Declare empty arrays to store specific values
radius = zeros(len(aperture))
signal = zeros(len(aperture))
background = zeros(len(aperture))
flux = zeros(len(aperture))

# Populate arrays with values
radius = aperture[:,8]
signal = aperture[:,6]
background = aperture[:,5]

for i in range(0, len(flux)):
    flux[i] = signal[i] - background[i]

# Plot Data
figure(1)
plot(radius, flux, 'bx')
xlabel('Radius of Aperture/Pixels')
ylabel('Background Subtracted Signal')
title('Signal as a Function of Aperture Radius to Determine \n Optimum Aperture Size')
savefig('aperture_qatar1b')