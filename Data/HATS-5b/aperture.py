'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Aperture Determination for HATS-5 b
'''

from numpy import *
from matplotlib.pyplot import *

# Read in data
aperture = genfromtxt('aperture.txt', dtype = 'float')
#aperture2 = genfromtxt('aperture2.txt', dtype = 'float')
#aperture3 = genfromtxt('aperture3.txt', dtype = 'float')

# Populate arrays with values
signal = aperture[where(aperture==1)[0],6]
background = aperture[where(aperture==2)[0],6]
radius = aperture[where(aperture==1)[0],8]
'''
signal2 = aperture2[where(aperture2==3)[0],6]
background2 = aperture2[where(aperture2==4)[0],6]
radius2 = aperture2[where(aperture2==3)[0],8]

signal3 = aperture3[where(aperture3==3)[0],6]
background3 = aperture3[where(aperture3==4)[0],6]
radius3 = aperture3[where(aperture3==3)[0],8]
'''
# Determine measured flux
flux = zeros(len(signal))
flux2 = zeros(len(signal))
flux3 = zeros(len(signal))

for i in range(0, len(flux)):
    flux[i] = signal[i] - background[i]
    #flux2[i] = signal2[i] - background2[i]
    #flux3[i] = signal3[i] - background3[i]

# Plot Data
figure(1)
plot(radius, flux, 'bx', label='Frame 2')
#plot(radius, flux2, 'rx', label = 'Frame 50')
#plot(radius, flux3, 'gx', label = 'Frame 110')
xlabel('Radius of Aperture/Pixels')
ylabel('Background Subtracted Signal')
title('Signal as a Function of Aperture Radius to Determine \n Optimum Aperture Size')
legend(loc='best')
savefig('aperture.png')