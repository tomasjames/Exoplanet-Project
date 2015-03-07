# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 13:08:38 2013

@author: tomasjames
Preliminary Project Preparation
Graphing for Prelim Investigations
"""

###############################################################################
#---------------------- Preliminary Project Work -----------------------------#
###############################################################################

#Import modules for code
import csv
from numpy import nan


############################### Read data #####################################

#Import data from exoplanet catalog csv file
with open('Catalogs/complete.csv', 'Ur') as openData:
    data = list(csv.reader(openData, delimiter=','))
 

############################# Section Data ####################################

#Import header
header = data[0]

#Declare empty lists to store different detection methods 
radial = []
imaging = []
transit = []
pulsar = []
astrometry = []
microlensing = []
    
#Sort data by detection method
for i in range(1, len(data)):
    if (data[i])[54] == 'Radial Velocity':
        radial.append(data[i])
  
    elif (data[i])[55] == 'Imaging':
        imaging.append(data[i])
  
    elif (data[i])[55] == 'Transit':
        transit.append(data[i])
        
    elif (data[i])[55] == 'pulsar':
        pulsar.append(data[i])
        
    elif (data[i])[55] == 'Astrometry':
        astrometry.append(data[i])
        
    elif (data[i])[55] == 'Microlensing':
        microlensing.append(data[i])
        
#Define arrays to store results
eccentricity = zeros(len(radial))     
mass = zeros(len(radial))
stellarMass = zeros(len(radial))
radius = zeros(len(radial))
period = zeros(len(radial))
effTemp = zeros(len(radial))
seperation = zeros(len(radial))

#Extract information into arrays
for i in range (1, len(radial)):

    if (data[1][13]) == ''

'''
    if (radial[i])[13] == '':
        eccentricity[i] = nan
    else:
        eccentricity[i-1] = float((radial[i])[13])
        
    if (radial[i])[1] == '':
        mass[i-1] = nan
    else:
        mass[i-1] = float((radial[i])[1])
        
    if (radial[i])[70] == '':
        stellarMass[i-1] = nan  
    else:
        stellarMass[i-1] = float((radial[i])[70])
        
    if (radial[i])[4] == '':
        radius[i] = nan
    else:
        radius[i-1] = float((radial[i])[4])
        
    if (radial[i])[7] == '':
        period[i-1] = nan
    else:
        period[i-1] = float((radial[i])[7])
        
    if (radial[i])[73] == '':
        effTemp[i-1] = nan
    else:
        effTemp[i-1] = float((radial[i])[73])
'''
########################## Initiate plotting #################################

'''        
#Plot graphs to investigate relations
figure(1)
subplot(2, 1, 1)
plot(eccentricity, mass, 'bx')
xlabel('eccentricity')
ylabel('mass')

subplot(2, 1, 2)
plot(log(eccentricity), log(mass), 'bx')
xlabel('eccentricity')
ylabel('mass')

figure(2)
plot(eccentricity, stellarMass, 'rx')
xlabel('eccentricity')
ylabel('stellarMass')

figure(3)
plot(effTemp, mass, 'gx')
xlabel('effTemp')
ylabel('mass')

figure(4)
subplot(2, 1, 1)
plot(stellarMass, mass, 'bx')
xlabel('stellarMass')
ylabel('mass')

subplot(2, 1, 2)
plot(log(stellarMass), log(mass), 'bx')
xlabel('stellarMass')
ylabel('mass')

figure(5)
hist(mass, 50)
xlabel('Mass')
ylabel('Frequency')
'''
