# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 13:48:31 2013

@author: tomasjames
"""

###############################################################################
#---------------------- Preliminary Project Work -----------------------------#
###############################################################################

#Import modules
import astropy
from astropy.table import Table, Column


############################### Read data #####################################

#Read votable using astropy
data = Table.read('Catalogs/complete.vot', format = 'votable')


############################# Section Data ####################################

#Declare empty lists to store different detection methods 


########################## Initiate plotting #################################

figure(1)
plot(transit[29], transit[1])


'''
figure(1)
subplot(2, 1, 1)
plot(data['col2'], data['col3'], 'bx')
xlabel('Effective Temperature of Host Star/K')
ylabel('Mass of Detected Exoplanet/MJup')

subplot(2, 1, 2)
plot(data['col30'], data['col2'], 'bx')
xlabel('Effective Temperature of Host Star/K')
ylabel('Mass of Detected Exoplanet/MJup')
xlim(2500, 7000)
ylim(0, 5)

suptitle('Mass of Detected Exoplanet as a Function of Host Star Eff. Temp.')

#Plot histograms for entire data set

figure(3)
subplot(2, 2, 1)
hist(data['col2'], bins = 70, range = (min(data['col2']), max(data['col2'])))
axvline(x = 1, color = 'g', linestyle = 'dashed', label = 'Jupiter Mass')
tight_layout()
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
title('Range: Full Data Set', fontsize = 12)
legend(loc='best')

subplot(2, 2, 2)
hist(data['col2'], bins = 60, range = (min(data['col2']), 10))
axvline(x = 1, color = 'g', linestyle = 'dashed', label = 'Jupiter Mass')
tight_layout()
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
title('Range: 0 - 10 MJup', fontsize = 12)
legend(loc='best')

subplot(2, 2, 3)
hist(data['col2'], bins = 50, range = (min(data['col2']), 1))
axvline(x = 0.003, color = 'g', linestyle = 'dashed', label = 'Earth Mass')
tight_layout()
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
title('Range: 0 - 1 MJup', fontsize = 12)
legend(loc='best')

subplot(2, 2, 4)
hist(data['col2'], bins = 50, range = (min(data['col2']), 0.1))
axvline(x = 0.003, color = 'g', linestyle = 'dashed', label = 'Earth Mass')
tight_layout()
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
title('Range: 0 - 0.1 MJup', fontsize = 12)
legend(loc='best')

suptitle('Histograms of Exoplanet Mass for Differing Mass Ranges for All Detections', fontsize = 14)

#Plot histograms for imaging data set

figure(4)
hist(imaging['col2'], bins = 20, range = (min(imaging['col2']), max(imaging['col2'])))
axvline(x = 1, color = 'g', linestyle = 'dashed', label = 'Jupiter Mass')
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
legend(loc='best')

title('Histograms of Exoplanet Mass for Imaging Detection', fontsize = 14)

#Plot histograms for microlensing data set

figure(5)
hist(microlensing['col2'], bins = 20, range = (min(microlensing['col2']), max(microlensing['col2'])))
axvline(x = 1, color = 'g', linestyle = 'dashed', label = 'Jupiter Mass')
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
legend(loc='best')

title('Histograms of Exoplanet Mass for Microlensing Detection', fontsize = 14)

#Plot histograms for radial data set

figure(6)
hist(radial['col2'], bins = 100, range = (min(radial['col2']), max(radial['col2'])))
axvline(x = 1, color = 'g', linestyle = 'dashed', label = 'Jupiter Mass')
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
legend(loc='best')

title('Histograms of Exoplanet Mass for Radial-Velocity Detection', fontsize = 14)

#Plot histograms for timing data set

figure(7)
hist(timing['col2'], bins = 30, range = (min(timing['col2']), max(timing['col2'])))
axvline(x = 1, color = 'g', linestyle = 'dashed', label = 'Jupiter Mass')
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
legend(loc='best')

title('Histograms of Exoplanet Mass for Timing Detection', fontsize = 14)

#Plot histograms for transiting data set

figure(8)
hist((transiting['col2']), bins = 70, range = (min(transiting['col2']), max(transiting['col2'])))
axvline(x = 1, color = 'g', linestyle = 'dashed', label = 'Jupiter Mass')
xlabel('Mass of Detected Exoplanet/Mjup')
ylabel('Frequency')
legend(loc='best')

title('Histograms of Exoplanet Mass for Transiting Detection', fontsize = 14)

#Plot graph of mass vs seperation
figure(9)
subplot(2, 1, 1)
tight_layout()
plot(data['col5'], data['col2'], 'bx')
xlabel('Seperation/AU')
ylabel('Mass of Detected Exoplanet/Mjup')
title('Range: Full Data Set', fontsize = 10)

subplot(2, 1, 2)
tight_layout()
plot(data['col5'], data['col2'], 'bx')
xlabel('Seperation/AU')
ylabel('Mass of Detected Exoplanet/Mjup')
xlim(0, 1)
ylim(0, 5)
title('Range: 0 - 1 Seperation, 0 - 5 MJup', fontsize = 10)

suptitle('Mass as a Function of Separation for Entire Catalog', fontsize = 14)

#Find total number of exoplanets at or below 1 MJup in mass
mass = array(data['col2'])

below = sum(i >= 1 for i in mass)

(below/len(mass))*100
'''