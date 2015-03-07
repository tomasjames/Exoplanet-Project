'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Test Transit Curve for WASP-12 b
'''

from os import *
from numpy import *
from matplotlib.pyplot import *

######### Index all files in directory containing photometry data #############

# Declare blank list to store filenames
fname = [] 

# Use walk function to index all filenames in given directory
for (dirpath, dirnames, filenames) in walk('Photometry'):
    fname.extend(filenames)
    break

# Above loop detects hidden files. Remove .DS_Store and aperture masks from list
fname.remove('.DS_Store')
fname.remove('GaiaPhotomIn.Dat')
fname.remove('GaiaPhotomOut.Dat')

# Determine quantities in .Dat files
test = genfromtxt('Photometry' + '/' + str(fname[0]), dtype='str')

# Declare 2d empty arrays to house imported values
# First row is star, second is background close to star
source = zeros([len(fname), 2]) 

calib1 = zeros([len(fname), 2])

calib2 = zeros([len(fname), 2])

calib3 = zeros([len(fname), 2])

# Loop through files in directory, take flux and background data into arrays
for i in range(0, len(fname)):
    f = genfromtxt('Photometry' + '/' + str(fname[i]), dtype='float')

    source[i, 0] = f[0,6]
    source[i, 1] = f[0,5]
    
    calib1[i, 0] = f[1,6]
    calib1[i, 1] = f[1,5]
    
    calib2[i, 0] = f[2,6]
    calib2[i, 1] = f[2,5]

    calib3[i, 0] = f[3,6]
    calib3[i, 1] = f[3,5]
    
# Calculate normalised flux from source
normalised1 = zeros(len(fname))
normalised2 = zeros(len(fname))
normalised3 = zeros(len(fname))

normalised1 = (source[:,0] - source[:,1])/(calib1[:,0] - calib1[:,1]) 
normalised2 = (source[:,0] - source[:,1])/(calib2[:,0] - calib2[:,1])  
normalised3 = (source[:,0] - source[:,1])/(calib3[:,0] - calib3[:,1])   

# Create array of time to plot against
time = linspace(1, len(fname), len(fname))

# Plot data
figure(1)
plot(time, normalised1, 'rx', label = 'Calibration 1')
xlabel('Time')  
ylabel('Normalised Flux')

legend(loc='best')

figure(2)
plot(time, normalised2, 'gx', label = 'Calibration 2')
xlabel('Time')  
ylabel('Normalised Flux')

legend(loc='best')

figure(3)
plot(time, normalised3, 'gx', label = 'Calibration 3')
xlabel('Time')  
ylabel('Normalised Flux')

legend(loc='best')