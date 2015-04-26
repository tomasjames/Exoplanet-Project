'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Transit Curve for HAT-P-25 b
'''

############################## Import modules #################################

from os import *
from numpy import *
from matplotlib.pyplot import *
from pyfits import *
from sklearn.metrics import mean_squared_error

# Import model
from modelv2 import *

############################### Read data #####################################

# Read photometry data
data = np.genfromtxt('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/HAT-P-25/Photometry/results.txt', dtype = 'float64')


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
calib5 = data[where(data==4)[0],6]
sky = data[where(data==7)[0],6]
radius = data[where(data==1)[0],8]


############################# Read FITS headers ###############################

# Declare list to house filenames
fname = []

# Walks through all files in /Raw/ and appends file name to fname if it 
# ends in .fits. Note: "../Raw/" steps up a directory to access /Raw/.
for file in listdir("/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/HAT-P-25/Raw"):
    if file.endswith(".fits"):
        fname.append('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/HAT-P-25/Raw/' + file)


########## Determine read noise and dark current from FITS headers ############

# Declare lists to house valies
rdnoise = 13.5 
dark = 0
pix = pi*(radius)**2

################################ Define constants for model ####################

Mstar = 1.01*1.9891E30 # kg
Rstar = 0.959*6.955E8 # m
Mplanet = 0.567*1.9E27 # kg
Rplanet = 1.19*(6.9E7) # m
radius = 4 # 4 pixel radii planet
end = 36537 # s
mid = 26257 # s
start = 16257 # s
nobs = 109 # Number of observations
dur = 36537-22157 # Frame 222 - Frame 40
a = 0.0466*149.6E9 # m
i = 86.0 # deg
mu = 0.522
app_mag = 11.69
sangle = 1./50. # Solid angle

'''
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
    error_source = sqrt((source+sky)+pix*(mean(dark)+mean(rdnoise)**2))
    error_calib = sqrt((calib+sky)+pix*(mean(dark)+mean(rdnoise)**2))
    error_sky = sqrt((sky)+pix*(mean(dark)+mean(rdnoise)**2)) 
    
    # Calculations the error in the curve equation
    error_curve = sqrt(((error_source**2)*(calib-source)**2 + (error_sky**2)*(source-calib)**2 + (error_calib**2)*(source-sky)**2)/(calib-sky)**4) 

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

delta_R = (20599516.011)/2

####################### Determine radius of exoplanet ########################

delta_F_1 = np.average(curve1[0:20]) - np.average(curve1[40:70])
delta_F_2 = np.average(curve2[0:20]) - np.average(curve2[40:70])
delta_F_3 = np.average(curve3[0:20]) - np.average(curve3[40:70])
delta_F_4 = np.average(curve4[0:20]) - np.average(curve4[40:70])
delta_F_5 = np.average(curve5[0:20]) - np.average(curve5[40:70])

det_Rplanet_1 = np.sqrt((Rstar**2)*(delta_F_1)*np.average(curve1[0:20]))
det_Rplanet_2 = np.sqrt((Rstar**2)*(delta_F_2))
det_Rplanet_3 = np.sqrt((Rstar**2)*(delta_F_3))
det_Rplanet_4 = np.sqrt((Rstar**2)*(delta_F_4))
det_Rplanet_5 = np.sqrt((Rstar**2)*(delta_F_5))

####################### Determine density of exoplanet ########################

rho_1 = (Mplanet/((4./3.)*np.pi*(det_Rplanet_1)**3))*(1000./(100)**3)
rho_2 = (Mplanet/((4./3.)*np.pi*(det_Rplanet_2)**3))*(1000./(100)**3)
rho_3 = (Mplanet/((4./3.)*np.pi*(det_Rplanet_3)**3))*(1000./(100)**3)
rho_4 = (Mplanet/((4./3.)*np.pi*(det_Rplanet_4)**3))*(1000./(100)**3)
rho_5 = (Mplanet/((4./3.)*np.pi*(det_Rplanet_5)**3))*(1000./(100)**3)

##################### Compute error in model to data fit #####################

error_fit_1 = np.sqrt(mean_squared_error(curve1, F))
error_fit_2 = np.sqrt(mean_squared_error(curve2, F))
error_fit_3 = np.sqrt(mean_squared_error(curve3, F))
error_fit_4 = np.sqrt(mean_squared_error(curve4, F))
error_fit_5 = np.sqrt(mean_squared_error(curve5, F))

##################### Compute error in calculated radius #####################

Rplanet_1_error = np.sqrt((delta_F_1/curve1[0])*(0.071)**2 + ((Rstar**2)/4)*(1/(delta_F_1/curve1[0]))*(0.005)**2 + ((Rstar**2)/4)*(delta_F_1/curve1[0]**3)*(error1[0])**2)
Rplanet_2_error = np.sqrt((delta_F_2/curve2[0])*(0.071)**2 + ((Rstar**2)/4)*(1/(delta_F_2/curve2[0]))*(0.005)**2 + ((Rstar**2)/4)*(delta_F_2/curve2[0]**3)*(error2[0])**2)
Rplanet_3_error = np.sqrt((delta_F_3/curve3[0])*(0.071)**2 + ((Rstar**2)/4)*(1/(delta_F_3/curve3[0]))*(0.005)**2 + ((Rstar**2)/4)*(delta_F_3/curve3[0]**3)*(error3[0])**2)
Rplanet_4_error = np.sqrt((delta_F_4/curve4[0])*(0.071)**2 + ((Rstar**2)/4)*(1/(delta_F_4/curve4[0]))*(0.005)**2 + ((Rstar**2)/4)*(delta_F_4/curve4[0]**3)*(error4[0])**2)
Rplanet_5_error = np.sqrt((delta_F_5/curve5[0])*(0.071)**2 + ((Rstar**2)/4)*(1/(delta_F_5/curve5[0]))*(0.005)**2 + ((Rstar**2)/4)*(delta_F_5/curve5[0]**3)*(error5[0])**2)

################ Determine error in density of exoplanet #####################

delta_rho_1 = np.sqrt((((3./(4*np.pi*(det_Rplanet_1)**3))**2)*((4.1756e25)**2)) + (((-9*Mplanet)/(4*np.pi*(det_Rplanet_1)**4))**2)*((Rplanet_1_error)**2))*(1000./(100)**3)
delta_rho_2 = np.sqrt((((3./(4*np.pi*(det_Rplanet_2)**3))**2)*((4.1756e25)**2)) + (((-9*Mplanet)/(4*np.pi*(det_Rplanet_2)**4))**2)*((Rplanet_2_error)**2))*(1000./(100)**3)
delta_rho_3 = np.sqrt((((3./(4*np.pi*(det_Rplanet_3)**3))**2)*((4.1756e25)**2)) + (((-9*Mplanet)/(4*np.pi*(det_Rplanet_3)**4))**2)*((Rplanet_3_error)**2))*(1000./(100)**3)
delta_rho_4 = np.sqrt((((3./(4*np.pi*(det_Rplanet_4)**3))**2)*((4.1756e26)**2)) + (((-9*Mplanet)/(4*np.pi*(det_Rplanet_4)**4))**2)*((Rplanet_4_error)**2))*(1000./(100)**3)
delta_rho_5 = np.sqrt((((3./(4*np.pi*(det_Rplanet_5)**3))**2)*((4.1756e25)**2)) + (((-9*Mplanet)/(4*np.pi*(det_Rplanet_5)**4))**2)*((Rplanet_5_error)**2))*(1000./(100)**3)

################################## Plot data ##################################

# Define array of frames to plot over
frames_data = linspace(1, len(source), len(source))
frames_model = linspace(1, 108, nobs)

figure(1)
#plot(frames, curve1, 'r.', label='Calibrated wrt calib1')
errorbar(frames_data, curve1, fmt='g.', yerr=error1, label='Calibrated wrt calib1')
plot(frames_model, F, 'b', label='EPTM')
#plot(frames, linspace(1,1,len(frames)))
ylim(0.96,1.02)
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for HAT-P-25')
text(20, 0.968, 'Radius: ('+str(det_Rplanet_1).format(1.0e9)+str(' +/- ')+str(Rplanet_1_error).format(1.0e9)+str(') m'))
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/HAT-P-25/Photometry/Graphs/curve1.png')



figure(2)
#plot(frames, curve2, 'r.', label='Calibrated wrt calib2')
errorbar(frames_data, curve2, fmt='', yerr=error2, label='Calibrated wrt calib2')
#plot(frames, linspace(1,1,len(frames)))
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for HAT-P-25')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/HAT-P-25/Photometry/Graphs/curve2.png')

figure(3)
#plot(frames_data, curve3, 'r.', label='Calibrated wrt calib3')
errorbar(frames_data, curve3, fmt='', yerr=error3)
#plot(frames, linspace(1,1,len(frames)))
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for HAT-P-25')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/HAT-P-25/Photometry/Graphs/curve3.png')

figure(4)
#plot(frames_data, curve4, 'r.', label='Calibrated wrt calib3')
errorbar(frames_data, curve4, fmt='', yerr=error4)
#plot(frames, linspace(1,1,len(frames)))
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for HAT-P-25')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/HAT-P-25/Photometry/Graphs/curve4.png')

figure(5)
#plot(frames_data, curve5, 'r.', label='Calibrated wrt calib3')
errorbar(frames_data, curve5, fmt='', yerr=error5)
#plot(frames, linspace(1,1,len(frames)))
xlabel('Frame Number')
ylabel('Calibrated Flux')
title('Transit Lightcurve for HAT-P-25')
legend(loc='best')
savefig('/Users/tomasjames/Documents/University/Cardiff/Project/Project/Data/HAT-P-25/Photometry/Graphs/curve5.png')

close('all')
