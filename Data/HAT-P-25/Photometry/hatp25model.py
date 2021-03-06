'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Model lightcurve for WASP-12b
'''

# Import model
from modelv2 import *

from matplotlib.pyplot import *

Mstar = 1.01*1.9891E30 # kg
Rstar = 0.959*6.955E8 # m
Mplanet = 0.567*1.9E27 # kg
Rplanet = 1.19*(6.9E7) # m
radius = 4 # 4 pixel radii planet
end = 36537 # s
mid = 26257 # s
start = 16257 # s
nobs = 240 # Number of observations
dur = 36537-22157 # Frame 222 - Frame 40
a = 0.0466*149.6E9 # m
i = 86.0 # deg
mu = 0.6
app_mag = 11.69
sangle = 1./42. # Solid angle

# Call model
F, X, T, R = model(start, mid, end, nobs, Mstar, Rstar, Mplanet, Rplanet, radius, a, mu, i, sangle)
#F, D = model(Rstar, a, mu, timestep, app_mag)

frames = np.linspace(1, 110, nobs)

# Plot data
figure(1)
plot(frames, F, label = 'Model for HAT-P-25')
xlabel('Distance from Transit Centre Point (m)')
ylabel('Relative Flux per Pixel Solid Angle')
ylim(0.970, 1.015)
legend(loc='best')
#savefig('Documents/University/Cardiff/Project/Project/Data/model_curve_hatp25.png', dpi=200)
show()

