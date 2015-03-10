'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Test Aperture Determination for HAT-P-25
'''

from numpy import *
from matplotlib.pyplot import *

# Read in data
aperture = genfromtxt('Documents/University/Cardiff/Project/Project/Data/HAT-P-25/aperture.txt', dtype = 'float')

# Declare empty arrays to store specific values
radius = zeros(len(aperture)/2)
signal = zeros(len(aperture)/2)
background = zeros(len(aperture)/2)
flux = zeros(len(aperture)/2)

# Populate arrays with values
signal = aperture[where(aperture==1)[0],6]
background = aperture[where(aperture==2)[0],6]
radius = aperture[where(aperture==1)[0],8]

for i in range(0, len(flux)):
    flux[i] = signal[i] - background[i]

# Plot Data
figure(1)
plot(radius, flux, 'bx')
xlabel('Radius of Aperture/Pixels')
ylabel('Background Subtracted Signal')
title('Signal as a Function of Aperture Radius to Determine \n Optimum Aperture Size')
savefig('aperture_hatp25b')
