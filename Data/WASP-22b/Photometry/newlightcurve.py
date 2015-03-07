'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Transit Curve for WASP-22 b
'''

############################## Import modules #################################

from os import *
from numpy import *
from matplotlib.pyplot import *
from pyfits import *


############################### Read data #####################################

# Read photometry data
data = genfromtxt('results.txt', dtype = 'float64')


##################### Extract quantities of interest ##########################

'''
The Numpy command where(data==x) finds the indices in the array data for which 
the element is equal to x. The appended [0] ensures that the result is 
1-dimensional - without this, the command returns the 2-dimensional position. 
The cumulative sum is found in column 6, which is the quantity of interest.

In this case, the index 1 indicates the target, source star whilst 7 indicates 
the sky aperture. Anything in between these two numbers is a calibration star.
'''

# Extract source, calibration star and sky in data
source = data[where(data==1)[0],6]
calib1 = data[where(data==2)[0],6]
calib2 = data[where(data==3)[0],6]
calib3 = data[where(data==4)[0],6]
calib4 = data[where(data==5)[0],6]
calib5 = data[where(data==6)[0],6]
sky = data[where(data==7)[0],6]
radius = data[where(data==1)[0],8]


############################# Read FITS headers ###############################

# Declare list to house filenames
fname = []

# Walks through all files in /Raw/ and appends file name to fname if it 
# ends in .fits. Note: "../Raw/" steps up a directory to access /Raw/.
for file in listdir("../Raw"):
    if file.endswith(".fits"):
        fname.append('../Raw/' + file)


########## Determine read noise and dark current from FITS headers ############

# Declare lists to house valies
rdnoise = []
dark = []
pix = pi*(radius)**2

# Loops over all values of filenames determined above and appends the value of
# the read noise RDNOISE and dark current DARKCURR to the lists
for i in range(0, len(fname)):
    rdnoise.append(getheader(fname[i])['RDNOISE'])
    dark.append(getheader(fname[i])['DARKCURR'])
    

########### Plot data to determine quality of calibration stars ###############

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
    ave = average(curve[0:10])
    
    # Normalises the curve using the baseline calculated earlier    
    norm = curve/ave

    # Determines the SNR of both the source and calib stars using the equation
    # for SNR from Handbook of CCD Astronomy
    error_source = sqrt((source+sky)+pix*(mean(dark)+mean(rdnoise)**2))
    error_calib = sqrt((calib+sky)+pix*(mean(dark)+mean(rdnoise)**2))
    error_sky = sqrt((sky)+pix*(mean(dark)+mean(rdnoise)**2)) 
    
    # Calculations the error in the curve equation
    error_curve = sqrt(((error_source**2)*(calib-source)**2 + (error_sky**2)*(source-calib)**2 + (error_calib**2)*(source-sky)**2)/(calib-sky)**4)
    
    return curve, error_curve
 
   
############### Run function to produce 3 calibrated data sets ################

curve1, error1 = lightcurve(source, calib1, sky, pix, dark, rdnoise)
curve2, error2 = lightcurve(source, calib2, sky, pix, dark, rdnoise)
curve3, error3 = lightcurve(source, calib3, sky, pix, dark, rdnoise)
curve4, error4 = lightcurve(source, calib4, sky, pix, dark, rdnoise)
curve5, error5 = lightcurve(source, calib5, sky, pix, dark, rdnoise)


################################## Plot data ##################################

# Define array of frames to plot over
frames = linspace(1, len(source), len(source))

figure(6)
errorbar(frames, curve1, fmt='', yerr=error1, label='Calibrated wrt calib1')
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('curve1.png')

figure(7)
errorbar(frames, curve2, fmt='', yerr=error2, label='Calibrated wrt calib2')
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('curve2.png')

figure(8)
errorbar(frames, curve3, fmt='', yerr=error3, label='Calibrated wrt calib3')
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('curve3.png')

figure(9)
errorbar(frames, curve4, fmt='', yerr=error4, label='Calibrated wrt calib4')
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('curve4.png')

figure(10)
errorbar(frames, curve5, fmt='', yerr=error5, label='Calibrated wrt calib5')
xlabel('Frame Number')
ylabel('Calibrated Flux')
legend(loc='best')
savefig('curve5.png')