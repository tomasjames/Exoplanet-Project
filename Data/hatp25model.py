'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Model lightcurve for WASP-12b
'''

# Import model
from modelv2 import *

from matplotlib.pyplot import *

Mstar = 1.01*1.9891e30 # kg
Rstar = 0.959*6.955e8 # m
Mplanet = 0.567*1.9e27 # kg
Rplanet = 1.19*(6.9e7) # m
radius = 5 # 4 pixel radii planet
end = 36537 # s
mid = 26257 # s
start = 16257 # s
nobs = 240 # Number of observations
dur = 36537-22157 # Frame 222 - Frame 40
a = 0.0466*149.6e9 # m
i = 86.0 # deg
mu = 1
app_mag = 11.69

# Call model
#x, y = model(Mstar, Rstar, Mplanet, Rplanet, radius, obs_end, trans_mid, nobs, a, i, mu, app_mag)
F, X, T = model(start, mid, end, nobs, Mstar, Rstar, Mplanet, Rplanet, radius, a, mu, 100)
#F, D = model(Rstar, a, mu, timestep, app_mag)

# Plot data
figure(1)
plot(X, F, label = 'Model for HAT-P-25')
xlabel('Distance from Transit Centre Point (m)')
ylabel('Incident Flux per Pixel Solid Angle (W/m**2/sr)')
legend(loc='best')
savefig('Documents/University/Cardiff/Project/Project/Data/model_curve_hatp25.png', dpi=200)
show()

