'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Model lightcurve for WASP-12b
'''

# Import model
from modelv2 import *

from matplotlib.pyplot import *

Mstar = 1.35*1.9891e30 # kg
Rstar = 1.599*6.955e8 # m
Mplanet = 1.404*1.9e27 # kg
Rplanet = 1.736*(6.9e7) # m
radius = 4 # 4 pixel radii planet
end = 36537 # s
mid = 26257 # s
start = 16257 # s
nobs = 240 # Number of observations
dur = 36537-22157 # Frame 222 - Frame 40
a = 0.02293*149.6e9 # m
i = 86.0 # deg
mu = 0.658
app_mag = 11.69

# Call model
#x, y = model(Mstar, Rstar, Mplanet, Rplanet, radius, obs_end, trans_mid, nobs, a, i, mu, app_mag)
F, X, D = model(start, mid, end, nobs, Mstar, Rstar, Mplanet, Rplanet, radius, a, mu)
#F, D = model(Rstar, a, mu, timestep, app_mag)

# Plot data
figure(1)
plot(X, F, label = 'Model for WASP-12b')
xlabel('Distance from Transit Centre Point (m)')
ylabel('Incident Flux per Pixel Solid Angle (W/m**2/sr)')
legend(loc='best')
savefig('Documents/University/Cardiff/Project/Project/Data/model_curve.png', dpi=200)
show()

