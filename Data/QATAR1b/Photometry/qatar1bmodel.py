'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Model lightcurve for WASP-12b
'''

# Import model
from modelv2 import *

from matplotlib.pyplot import *

Mstar = 0.85*1.9891e30 # kg
Rstar = 0.823*6.955e8 # m
Mplanet = 1.09*1.9e27 # kg
Rplanet = 1.164*(6.9e7) # m
radius = 4 # 4 pixel radii planet
end = 29108 # s
mid = 22724 # s
start = 16340 # s
nobs = 181 # Number of observations
a = 0.02343*149.6e9 # m
i = 86.0 # deg
mu = 0.6
sangle = 1./36.

# Call model
#x, y = model(Mstar, Rstar, Mplanet, Rplanet, radius, obs_end, trans_mid, nobs, a, i, mu, app_mag)
F, X, F_A, R = model(start, mid, end, nobs, Mstar, Rstar, Mplanet, Rplanet, radius, a, mu, 100, sangle)
#F, D = model(Rstar, a, mu, timestep, app_mag)

# Plot data
figure(1)
plot(X, F, label = 'Model for QATAR-1b')
xlabel('Distance from Transit Centre Point (m)')
ylabel('Incident Flux per Pixel Solid Angle (W/m**2/sr)')
ylim(0.95, 1.03)
legend(loc='best')
savefig('Documents/University/Cardiff/Project/Project/Data/model_curve_qatar1b.png', dpi=200)
show()

