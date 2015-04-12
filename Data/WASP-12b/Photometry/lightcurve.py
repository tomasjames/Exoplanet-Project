'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Transit Curve for WASP-12 b
'''

############################## Import modules #################################

from os import *
import numpy as np
from matplotlib.pyplot import *
from pyfits import *
from modelv2 import *


############################### Read data #####################################

# Read photometry data
data = np.genfromtxt('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-12b/Photometry/results.txt', dtype = 'float64')


##################### Extract quantities of interest ##########################

'''
The Numpy command where(data==x) finds the indices in the array data for which 
the element is equal to x. The appended [0] ensures that the result is 
1-dimensional - without this, the command returns the 2-dimensional position. 
The cumulative sum is found in column 6, which is the quantity of interest.

In this case, the index 1 indicates the target, source star whilst 5 indicates 
the sky aperture. Anything in between these two numbers is a calibration star.
'''

# Extract source, calibration star and sky in data
source = data[np.where(data==1)[0],6]
calib1 = data[np.where(data==2)[0],6]
calib2 = data[np.where(data==3)[0],6]
calib3 = data[np.where(data==4)[0],6]
calib4 = data[np.where(data==5)[0],6]
calib5 = data[np.where(data==6)[0],6]
sky = data[np.where(data==7)[0],6]
radius = data[np.where(data==1)[0],8]


############################# Read FITS headers ###############################

# Declare list to house filenames
fname = []

# Walks through all files in /Raw/ and appends file name to fname if it 
# ends in .fits. Note: "../Raw/" steps up a directory to access /Raw/.
for file in listdir("/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-12b/Raw"):
    if file.endswith(".fits"):
        fname.append('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-12b/Raw/' + file)


########## Determine read noise and dark current from FITS headers ############

# Declare lists to house valies
rdnoise = 13.5 
dark = 0
pix = np.pi*(radius)**2

############################# Declare constants ################################

Mstar = 1.35*1.9891e30 # kg
Rstar = 1.599*6.955e8 # m
Mplanet = 1.404*1.9e27 # kg
Rplanet = 1.736*(6.9e7) # m
radius = 4. # 4 pixel radii planet
end = 35200 # s
mid = 26257 # s
start = 17600 # s
nobs = 240 # Number of observations
a = 0.02293*149.6e9 # m
i = 86.0 # deg
mu = 0.49
app_mag = 11.69
sangle = 1./65.

########### Plot data to determine quality of calibration stars ###############
'''
figure(1)

subplot(2, 2, 1)
plot((calib1-sky)/(calib2-sky), label='(calib1-sky)/(calib2-sky)')
ylabel('Background Subtracted Count')
legend(loc='best')

subplot(2, 2, 2)
plot((calib1-sky)/(calib3-sky), label='(calib1-sky)/(calib3-sky)')
ylabel('Background Subtracted Count')
legend(loc='best')

subplot(2, 2, 3)
plot((calib1-sky)/(calib4-sky), label='(calib1-sky)/(calib4-sky)')
xlabel('Frame Number')
ylabel('Background Subtracted Count')
legend(loc='best')

subplot(2, 2, 4)
plot((calib1-sky)/(calib5-sky), label='(calib1-sky)/(calib5-sky)')
xlabel('Frame Number')
ylabel('Background Subtracted Count')
legend(loc='best')

savefig('calibration1.png')


figure(2)

subplot(2, 2, 1)
plot((calib2-sky)/(calib3-sky), label='(calib2-sky)/(calib3-sky)')
ylabel('Background Subtracted Count')
legend(loc='best')

subplot(2, 2, 2)
plot((calib2-sky)/(calib4-sky), label='(calib2-sky)/(calib4-sky)')
ylabel('Background Subtracted Count')
legend(loc='best')

subplot(2, 2, 3)
plot((calib2-sky)/(calib5-sky), label='(calib2-sky)/(calib5-sky)')
ylabel('Background Subtracted Count')
legend(loc='best')

savefig('calibration2.png')

figure(3)

subplot(2, 1, 1)
plot((calib3-sky)/(calib4-sky), label='(calib3-sky)/(calib4-sky)')
ylabel('Background Subtracted Count')
legend(loc='best')

subplot(2, 1, 2)
plot((calib3-sky)/(calib5-sky), label='(calib3-sky)/(calib5-sky)')
ylabel('Background Subtracted Count')
legend(loc='best')

savefig('calibration3.png')

figure(4)

plot((calib4-sky)/(calib5-sky), label='(calib4-sky)/(calib5-sky)')
ylabel('Background Subtracted Count')
legend(loc='best')
'''

########################## Define lightcurve equation #########################

def lightcurve(source, calib, sky, pix, dark, rdnoise):
   
    '''
    Uses the simple lightcurve equation to calibrate a source star in order to 
    observe a transit (if one is present). Takes data for the source star, one
    calibration star and a sky data. 
    
    Also normalises the transit by dividing all points by an off-transit 
    average.
    '''
    
    # Calculates lightcurve values 
    curve = (source-sky)/(calib-sky)

    # Takes the first and last 10 off-transit values and averages 
    # them for a baseline    
    ave = (np.average(curve[0:20]) + np.average(curve[-20:-1]))/2
    
    # Normalises the curve using the baseline calculated earlier    
    norm = curve/ave

    # Determines the SNR of both the source and calib stars using the equation
    # for SNR from Handbook of CCD Astronomy
    error_source = np.sqrt((source+sky)+pix*(np.mean(dark)+np.mean(rdnoise)**2))
    error_calib = np.sqrt((calib+sky)+pix*(np.mean(dark)+np.mean(rdnoise)**2))
    error_sky = np.sqrt((sky)+pix*(np.mean(dark)+np.mean(rdnoise)**2)) 
    
    # Calculations the error in the curve equation
    error_curve = np.sqrt(((error_source**2)*(calib-source)**2 + (error_sky**2)*(source-calib)**2 + (error_calib**2)*(source-sky)**2)/(calib-sky)**4) 

    error_norm  = error_curve/ave
    
    #return curve, error_curve
    return norm, error_norm
    
############### Run function to produce 3 calibrated data sets ################

curve1, error1 = lightcurve(source, calib1, sky, pix, dark, rdnoise)
curve2, error2 = lightcurve(source, calib2, sky, pix, dark, rdnoise)
curve3, error3 = lightcurve(source, calib3, sky, pix, dark, rdnoise)
curve4, error4 = lightcurve(source, calib4, sky, pix, dark, rdnoise)
curve5, error5 = lightcurve(source, calib5, sky, pix, dark, rdnoise)

################################## Call model #################################

F, X, T, R = model(start, mid, end, nobs, Mstar, Rstar, Mplanet, Rplanet, radius, a, mu, i, sangle)

################################## Plot data ##################################

# Define array of frames to plot over
frames = np.linspace(1, len(source), len(source))
frames_model = np.linspace(1, 220, nobs)

figure(6)
#plot(frames, curve1, 'r.', label='Calibrated wrt calib1')
errorbar(frames, curve1, fmt='', yerr=error1, label='Calibrated wrt calib1')
#plot(frames, linspace(1,1,len(frames)))
plot(frames_model, F)
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-12b/Photometry/Graphs/curve1.png')

det_Rplanet_1 = np.sqrt((Rstar**2)*(np.average(curve1[0:20]) - min(curve1)))

figure(7)
#plot(frames, curve2, 'r.', label='Calibrated wrt calib2')
errorbar(frames, curve2, fmt='', yerr=error2, label='Calibrated wrt calib2')
plot(frames_model, F)
#plot(frames, linspace(1,1,len(frames)))
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-12b/Photometry/Graphs/curve2.png')

det_Rplanet_2 = np.sqrt((Rstar**2)*(np.average(curve2[0:20]) - min(curve2)))

figure(8)
#plot(frames, curve3, 'r.', label='Calibrated wrt calib3')
errorbar(frames, curve3, fmt='', yerr=error3, label='Calibrated wrt calib3')
plot(frames_model, F)
#plot(frames, linspace(1,1,len(frames)))
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-12b/Photometry/Graphs/curve3.png')

det_Rplanet_3 = np.sqrt((Rstar**2)*(np.average(curve3[0:20]) - min(curve3)))

figure(9)
#plot(frames, curve4, 'r.', label='Calibrated wrt calib4')
errorbar(frames, curve4, fmt='', yerr=error4, label='Calibrated wrt calib4')
plot(frames_model, F)
#plot(frames, linspace(1,1,len(frames)))
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-12b/Photometry/Graphs/curve4.png')

det_Rplanet_4 = np.sqrt((Rstar**2)*(np.average(curve4[0:20]) - min(curve4)))

figure(10)
#plot(frames, curve5, 'r.', label='Calibrated wrt calib5')
errorbar(frames, curve5, fmt='', yerr=error5, label='Calibrated wrt calib5')
plot(frames_model, F)
#plot(frames, linspace(1,1,len(frames)))
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-12b/Photometry/Graphs/curve5.png')

det_Rplanet_5 = np.sqrt((Rstar**2)*(np.average(curve5[0:20]) - min(curve5)))

close('all')
