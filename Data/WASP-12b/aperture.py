'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Aperture Determination for WASP-12 b
'''

from numpy import *
from matplotlib.pyplot import *

# Read in data
aperture = genfromtxt('aperture.txt', dtype = 'float')

# Populate arrays with values
signal = aperture[where(aperture==1)[0],6]
background = aperture[where(aperture==2)[0],6]
radius = aperture[where(aperture==1)[0],8]

# Determine measured flux
flux = zeros(len(signal))

for i in range(0, len(flux)):
    flux[i] = signal[i] - background[i]

# Plot Data
figure(1)
plot(radius, flux, 'bx', label='Frame 1')
xlabel('Radius of Aperture/Pixels')
ylabel('Background Subtracted Signal')
title('Signal as a Function of Aperture Radius to Determine \n Optimum Aperture Size')
legend(loc='best')
savefig('aperture.png')