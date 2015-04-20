'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Transit Curve for WASP-78 b
'''

############################## Import modules #################################

from os import *
from numpy import *
from matplotlib.pyplot import *

############################### Read data #####################################

# Read photometry data
data = genfromtxt('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/results.txt', dtype = 'float64')


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
source = data[where(data==1)[0],6]
calib1 = data[where(data==2)[0],6]
calib2 = data[where(data==3)[0],6]
calib3 = data[where(data==4)[0],6]
calib4 = data[where(data==5)[0],6]
calib5 = data[where(data==6)[0],6]
calib6 = data[where(data==7)[0],6]
calib7 = data[where(data==8)[0],6]
calib8 = data[where(data==9)[0],6]
calib9 = data[where(data==10)[0],6]
calib10 = data[where(data==11)[0],6]
calib11 = data[where(data==12)[0],6]
calib12 = data[where(data==13)[0],6]
calib13 = data[where(data==14)[0],6]
sky = data[where(data==15)[0],6]
radius = data[where(data==1)[0],8]

############################# Read FITS headers ###############################

# Declare list to house filenames
fname = []
'''
# Walks through all files in /Raw/ and appends file name to fname if it 
# ends in .fits. Note: "../Raw/" steps up a directory to access /Raw/.
for file in listdir("../Raw"):
    if file.endswith(".fits"):
        fname.append('../Raw/' + file)
'''

########## Determine read noise and dark current from FITS headers ############

# Declare lists to house valies
rdnoise = 13.5 
dark = 0
pix = pi*(radius)**2

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

savefig('calibration4.png')


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
curve6, error6 = lightcurve(source, calib6, sky, pix, dark, rdnoise)
curve7, error7 = lightcurve(source, calib7, sky, pix, dark, rdnoise)
curve8, error8 = lightcurve(source, calib8, sky, pix, dark, rdnoise)
curve9, error9 = lightcurve(source, calib9, sky, pix, dark, rdnoise)
curve10, error10 = lightcurve(source, calib10, sky, pix, dark, rdnoise)
curve11, error11 = lightcurve(source, calib11, sky, pix, dark, rdnoise)
curve12, error12 = lightcurve(source, calib12, sky, pix, dark, rdnoise)
curve13, error13 = lightcurve(source, calib13, sky, pix, dark, rdnoise)


################################## Plot data ##################################

# Define array of frames to plot over
frames = linspace(1, len(source), len(source))

figure(5)
errorbar(frames, curve1, yerr=error1)
plot(frames, curve1, 'r.', label='Calibrated wrt calib1')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve1.png')

figure(6)
errorbar(frames, curve2, yerr=error2)
plot(frames, curve2, 'r.', label='Calibrated wrt calib2')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve2.png')

figure(7)
errorbar(frames, curve3, yerr=error3)
plot(frames, curve3, 'r.', label='Calibrated wrt calib3')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve3.png')

figure(8)
errorbar(frames, curve4, yerr=error4)
plot(frames, curve4, 'r.', label='Calibrated wrt calib4')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve4.png')

figure(9)
errorbar(frames, curve5, yerr=error5)
plot(frames, curve5, 'r.', label='Calibrated wrt calib5')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve5.png')

figure(10)
errorbar(frames, curve6, yerr=error6)
plot(frames, curve6, 'r.', label='Calibrated wrt calib6')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve6.png')

figure(11)
errorbar(frames, curve7, yerr=error7)
plot(frames, curve7, 'r.', label='Calibrated wrt calib7')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve7.png')

figure(12)
errorbar(frames, curve8, yerr=error8)
plot(frames, curve8, 'r.', label='Calibrated wrt calib8')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve8.png')

figure(13)
errorbar(frames, curve9, yerr=error9)
plot(frames, curve9, 'r.', label='Calibrated wrt calib9')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve9.png')

figure(14)
errorbar(frames, curve10, yerr=error10)
plot(frames, curve10, 'r.', label='Calibrated wrt calib10')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve10.png')

figure(15)
errorbar(frames, curve11, yerr=error11)
plot(frames, curve11, 'r.', label='Calibrated wrt calib11')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve11.png')

figure(16)
errorbar(frames, curve12, yerr=error12)
plot(frames, curve12, 'r.', label='Calibrated wrt calib12')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve12.png')

figure(17)
errorbar(frames, curve13, yerr=error13)
plot(frames, curve13, 'r.', label='Calibrated wrt calib13')
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for WASP-78b')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/WASP-78b/Photometry/curve13.png')
