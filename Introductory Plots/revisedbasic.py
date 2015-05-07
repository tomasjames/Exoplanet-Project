# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 16:41:44 2014

@author: tomasjames
"""

###############################################################################
#---------------------- Preliminary Project Work -----------------------------#
###############################################################################

#Import modules
from astropy.table import Table, Column
from matplotlib.pyplot import *
from numpy import *

############################### Read data #####################################

#Import data
data = Table.read('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Introductory Plots/Catalogs/complete.xml',format='votable')

###################### Extract variables from data ############################

#Declare empty lists to store data
radial_mass = []
radial_radius = []
radial_major = []
radial_smass = []
radial_period = []
radial_inc = []

transit_mass = []
transit_radius = []
transit_major = []
transit_smass = []
transit_period = []
transit_inc = []

imaging_mass = []
imaging_radius = []
imaging_major = []
imaging_smass = []
imaging_period = []
imaging_inc = []

astrometry_mass = []
astrometry_radius = []
astrometry_major = []
astrometry_smass = []
astrometry_period = []
astrometry_inc = []

microlensing_mass = []
microlensing_radius = []
microlensing_major = []
microlensing_smass = []
microlensing_period = []
microlensing_inc = []

pulsar_mass = []
pulsar_radius = []
pulsar_major = []
pulsar_smass = []
pulsar_period = []
pulsar_inc = []

ttv_mass = []
ttv_radius = []
ttv_major = []
ttv_smass = []
ttv_period = []
ttv_inc = []

#Iterate through detection type column and section into lists
for i in range(0, len(data['col55'])):
    if (data['col55'])[i] == 'Radial Velocity':
        radial_mass.append((data['col2'])[i])
        radial_radius.append((data['col5'])[i])
        radial_major.append((data['col11'])[i])
        radial_smass.append((data['col70'])[i])
        radial_period.append((data['col8'])[i])
        radial_inc.append((data['col17'])[i])        

    elif (data['col55'])[i] == 'Transit':
        transit_mass.append((data['col2'])[i])
        transit_radius.append((data['col5'])[i])
        transit_major.append((data['col11'])[i])
        transit_smass.append((data['col70'])[i])
        transit_period.append((data['col8'])[i])
        transit_inc.append((data['col17'])[i])        
        
    elif (data['col55'])[i] == 'Imaging':
        imaging_mass.append((data['col2'])[i])
        imaging_radius.append((data['col5'])[i])
        imaging_major.append((data['col11'])[i])
        imaging_smass.append((data['col70'])[i])        
        imaging_period.append((data['col8'])[i])
        imaging_inc.append((data['col17'])[i])

    elif (data['col55'])[i] == 'Astrometry':
        astrometry_mass.append((data['col2'])[i])
        astrometry_radius.append((data['col5'])[i])
        astrometry_major.append((data['col11'])[i])
        astrometry_smass.append((data['col70'])[i])
        astrometry_period.append((data['col8'])[i])
        astrometry_inc.append((data['col17'])[i])

    elif (data['col55'])[i] == 'Microlensing':
        microlensing_mass.append((data['col2'])[i])
        microlensing_radius.append((data['col5'])[i])
        microlensing_major.append((data['col11'])[i])
        microlensing_smass.append((data['col70'])[i])
        microlensing_period.append((data['col8'])[i])
        microlensing_inc.append((data['col17'])[i])
        
    elif (data['col55'])[i] == 'Pulsar':
        pulsar_mass.append((data['col2'])[i])
        pulsar_radius.append((data['col5'])[i])
        pulsar_major.append((data['col11'])[i])
        pulsar_smass.append((data['col70'])[i])
        pulsar_period.append((data['col8'])[i])
        pulsar_inc.append((data['col17'])[i])
        
    elif (data['col55'])[i] == 'TTV':
        ttv_mass.append((data['col2'])[i])
        ttv_radius.append((data['col5'])[i])
        ttv_major.append((data['col11'])[i])
        ttv_smass.append((data['col70'])[i])
        ttv_period.append((data['col8'])[i])
        ttv_inc.append((data['col17'])[i])

#Detect elements with empty cells and replace with nan to avoid plot errors      
for i in range(0, len(radial_mass)):
    
    if radial_mass[i] == '0':
        radial_mass[i] = nan
    
    if radial_radius[i] == '0':
        radial_radius[i] = nan
        
    if radial_major[i] == '0':
        radial_major[i] = nan
        
    if radial_smass[i] == '0':
        radial_smass[i] = nan
        
    if radial_period[i] == '0':
        radial_period[i] = nan
        
    if radial_inc[i] == '0':
        radial_inc[i] = nan
        
#Convert lists to arrays in order to plot histograms
transit_mass_arr = asarray(transit_mass)

#Remove nan elements
transit_mass_arr = transit_mass_arr[~isnan(transit_mass_arr)]

#Determine smallest and biggest masses
transit_mass_min = min(transit_mass_arr)
transit_mass_max = max(transit_mass_arr)

######################## Define Solar System data #############################

jupiter_radius = 1
earth_radius = 6371./69911

jupiter_mass = 1
earth_mass = 6.e24/1.9e27
sun_mass = 1

earth_major = 1.00000011
jupiter_major = 5.20336301

earth_period = 365
jupiter_period = 4332
    
############################### Plot data #####################################

earthline = ones(len(radial_mass))*earth_radius
line1 = linspace(1, 30, len(radial_mass))

#Open figure
figure(1)

#Plot mass-radius relation for all detection types
plot(radial_radius, radial_mass, 'm.', label = 'Radial Velocity: ' +str(len(radial_radius)))
plot(transit_radius, transit_mass, 'b.', label = 'Transit: ' +str(len(transit_radius)))
plot(imaging_radius, imaging_mass, 'r.', label = 'Imaging: ' +str(len(imaging_radius)))
plot(astrometry_radius, astrometry_mass, 'y.', label = 'Astrometry: ' +str(len(astrometry_radius)))
plot(microlensing_radius, microlensing_mass, 'g.', label = 'Microlensing: ' +str(len(microlensing_radius)))
plot(pulsar_radius, pulsar_mass, 'c.', label = 'Pulsar: ' +str(len(pulsar_radius)))
plot(ttv_radius, ttv_mass, 'k.', label = 'TTV: ' +str(len(ttv_radius)))
plot(earth_radius, earth_mass, 'g<', markersize = 10, label = 'Earth')
plot(jupiter_radius, jupiter_mass, 'r>', markersize = 10, label = 'Jupiter')
plot(earthline, line1, 'g--')

xlabel('$Exoplanet Radius /R_{Jup}$')
ylabel('$Exoplanetary Mass /M_{Jup}$')
title('Plot of Exoplanetary Mass as a Function of Radius')
legend(loc='best', prop={'size':10})

savefig('radius_mass')


line2 = linspace(10**(-3), 10**(2), len(radial_mass))

figure(2)
loglog(radial_radius, radial_mass, 'm.', label = 'Radial Velocity: ' +str(len(radial_radius)))
loglog(transit_radius, transit_mass, 'b.', label = 'Transit: ' +str(len(transit_radius)))
loglog(imaging_radius, imaging_mass, 'r.', label = 'Imaging: ' +str(len(imaging_radius)))
loglog(astrometry_radius, astrometry_mass, 'y.', label = 'Astrometry: ' +str(len(astrometry_radius)))
loglog(microlensing_radius, microlensing_mass, 'g.', label = 'Microlensing: ' +str(len(microlensing_radius)))
loglog(pulsar_radius, pulsar_mass, 'c.', label = 'Pulsar: ' +str(len(pulsar_radius)))
loglog(ttv_radius, ttv_mass, 'k.', label = 'TTV: ' +str(len(ttv_radius)))
loglog(earth_radius, earth_mass, 'g<', markersize = 10, label = 'Earth')
loglog(jupiter_radius, jupiter_mass, 'r>', markersize = 10, label = 'Jupiter')
plot(log(earthline), line2, 'g--')

xlabel('$Log(Exoplanet Radius) /R_{Jup}$')
ylabel('$Log(Exoplanetary Mass) /M_{Jup}$')
title('Logarithmic Plot of Exoplanetary Mass as a Function of Radius')
legend(loc='best', prop={'size':10})

savefig('log_radius_mass')


figure(3)
plot(radial_major, radial_smass, 'm.', label = 'Radial Velocity: ' +str(len(radial_major)))
plot(transit_major, transit_smass, 'b.', label = 'Transit: ' +str(len(transit_major)))
plot(imaging_major, imaging_smass, 'r.', label = 'Imaging: ' +str(len(imaging_major)))
plot(astrometry_major, astrometry_smass, 'y.', label = 'Astrometry: ' +str(len(astrometry_major)))
plot(microlensing_major, microlensing_smass, 'g.', label = 'Microlensing: ' +str(len(microlensing_major)))
plot(pulsar_major, pulsar_smass, 'c.', label = 'Pulsar: ' +str(len(pulsar_major)))
plot(ttv_major, ttv_smass, 'k.', label = 'TTV: ' +str(len(ttv_major)))
plot(earth_major, sun_mass, 'g<', markersize = 10, label = 'Earth')
plot(jupiter_major, sun_mass, 'r>', markersize = 10, label = 'Jupiter')

xlabel('$Semi-Major Axis / AU $')
ylabel('$Mass of Host Star / M_{Sun}$')
title('Mass of Host Star as a Function of Exoplanetary Semi-Major Axis')
legend(loc='best')

savefig('major_smass')


figure(4)
loglog(radial_major, radial_smass, 'm.', label = 'Radial Velocity: ' +str(len(radial_major)))
loglog(transit_major, transit_smass, 'b.', label = 'Transit: ' +str(len(transit_major)))
loglog(imaging_major, imaging_smass, 'r.', label = 'Imaging: ' +str(len(imaging_major)))
loglog(astrometry_major, astrometry_smass, 'y.', label = 'Astrometry: ' +str(len(astrometry_major)))
loglog(microlensing_major, microlensing_smass, 'g.', label = 'Microlensing: ' +str(len(microlensing_major)))
loglog(pulsar_major, pulsar_smass, 'c.', label = 'Pulsar: ' +str(len(pulsar_major)))
loglog(ttv_major, ttv_smass, 'k.', label = 'TTV: ' +str(len(ttv_major)))
loglog(earth_major, sun_mass, 'g<', markersize = 10, label = 'Earth')
loglog(jupiter_major, sun_mass, 'r>', markersize = 10, label = 'Jupiter')

xlabel('$Log(Semi-Major Axis) / AU $')
ylabel('$Log(Mass of Host Star) / M_{Sun}$')
title('Logarithmic Plot of Mass of Host Star as \n a Function of Exoplanetary Semi-Major Axis')
legend(loc='best')

savefig('log_major_smass')


figure(5)
plot(radial_period, radial_smass, 'm.', label = 'Radial Velocity: ' +str(len(radial_smass)))
plot(transit_period, transit_smass, 'b.', label = 'Transit: ' +str(len(transit_smass)))
plot(imaging_period, imaging_smass, 'r.', label = 'Imaging: ' +str(len(imaging_smass)))
plot(astrometry_period, astrometry_smass, 'y.', label = 'Astrometry: ' +str(len(astrometry_smass)))
plot(microlensing_period, microlensing_smass, 'g.', label = 'Microlensing: ' +str(len(microlensing_smass)))
plot(pulsar_period, pulsar_smass, 'c.', label = 'Pulsar: ' +str(len(pulsar_smass)))
plot(ttv_period, ttv_smass, 'k.', label = 'TTV: ' +str(len(ttv_smass)))
plot(earth_period, sun_mass, 'g<', markersize = 10, label = 'Earth')
plot(jupiter_period, sun_mass, 'r>', markersize = 10, label = 'Jupiter')

xlabel('$Exoplanet Period / Days $')
ylabel('$Mass of Host Star / M_{Sun}$')
title('Mass of Host Star as a Function of Exoplanet Period')
legend(loc='best')

savefig('period_smass')


figure(6)
loglog(radial_period, radial_smass, 'm.', label = 'Radial Velocity: ' +str(len(radial_smass)))
loglog(transit_period, transit_smass, 'b.', label = 'Transit: ' +str(len(transit_smass)))
loglog(imaging_period, imaging_smass, 'r.', label = 'Imaging: ' +str(len(imaging_smass)))
loglog(astrometry_period, astrometry_smass, 'y.', label = 'Astrometry: ' +str(len(astrometry_smass)))
loglog(microlensing_period, microlensing_smass, 'g.', label = 'Microlensing: ' +str(len(microlensing_smass)))
loglog(pulsar_period, pulsar_smass, 'c.', label = 'Pulsar: ' +str(len(pulsar_smass)))
loglog(ttv_period, ttv_smass, 'k.', label = 'TTV: ' +str(len(ttv_smass)))
loglog(earth_period, sun_mass, 'g<', markersize = 10, label = 'Earth')
loglog(jupiter_period, sun_mass, 'r>', markersize = 10, label = 'Jupiter')

xlabel('$Log(Period) / Days $')
ylabel('$Log(Mass of Host Star) / M_{Sun}$')
title('Logarithmic Plot of Mass of Host Star as a Function of Exoplanet Period')
legend(loc='best')

savefig('log_period_smass')


figure(7)
plot(radial_major, radial_period, 'm.', label = 'Radial Velocity: ' +str(len(radial_period)))
plot(transit_major, transit_period, 'b.', label = 'Transit: ' +str(len(transit_period)))
plot(imaging_major, imaging_period, 'r.', label = 'Imaging: ' +str(len(imaging_period)))
plot(astrometry_major, astrometry_period, 'y.', label = 'Astrometry: ' +str(len(astrometry_period)))
plot(microlensing_major, microlensing_period, 'g.', label = 'Microlensing: ' +str(len(microlensing_period)))
plot(pulsar_major, pulsar_period, 'c.', label = 'Pulsar: ' +str(len(pulsar_period)))
plot(ttv_major, ttv_period, 'k.', label = 'TTV: ' +str(len(ttv_period)))
plot(earth_major, earth_period, 'g<', markersize = 10, label = 'Earth')
plot(jupiter_major, jupiter_period, 'r>', markersize = 10, label = 'Jupiter')

xlabel('$Semi-Major Axis / AU$')
ylabel('$Exoplanet Period / Days $')
title('Mass of Host Star as a Function of Exoplanet Period')
legend(loc='best', prop={'size':10})

savefig('period_major')


figure(8)
loglog(radial_major, radial_period, 'm.', label = 'Radial Velocity: ' +str(len(radial_period)))
loglog(transit_major, transit_period, 'b.', label = 'Transit: ' +str(len(transit_period)))
loglog(imaging_major, imaging_period, 'r.', label = 'Imaging: ' +str(len(imaging_period)))
loglog(astrometry_major, astrometry_period, 'y.', label = 'Astrometry: ' +str(len(astrometry_period)))
loglog(microlensing_major, microlensing_period, 'g.', label = 'Microlensing: ' +str(len(microlensing_period)))
loglog(pulsar_major, pulsar_period, 'c.', label = 'Pulsar: ' +str(len(pulsar_period)))
loglog(ttv_major, ttv_period, 'k.', label = 'TTV: ' +str(len(ttv_period)))
loglog(earth_major, earth_period, 'g<', markersize = 10, label = 'Earth')
loglog(jupiter_major, jupiter_period, 'r>', markersize = 10, label = 'Jupiter')

xlabel('$Log(Semi-Major Axis) / AU$')
ylabel('$Log(Exoplanet Period) / Days $')
title('Logarithmic Plot of Mass of Host Star as a Function of Exoplanet Period')
legend(loc='best', prop={'size':10})

savefig('log_period_major')

figure(9)
hist(transit_mass_arr, bins=40)

xlabel('$Exoplanetary Mass /M_{Jup}$')
ylabel('$Frequency$')
title('Histogram Showing the Distribution of Exoplanetary Masses \n Detected Using Transit Method')

savefig('hist_transit_mass')
