'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Lightcurve production script
'''

from numpy import *

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

    error_norm  = error_curve/ave
    
    #return curve, error_curve
    return norm, error_norm
