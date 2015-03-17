'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Prototype for modelling exoplanetary lightcurves
'''

import numpy as np
from matplotlib.pyplot import *
from itertools import product

def model(start, mid, end, nobs, Mstar, Rstar, Mplanet, Rplanet, radius, a, mu, npix):

    '''
    Takes a planetary system consisting of an exoplanet of mass M_planet
    orbiting a star of mass M_star and radius R_star and models a light-
    curve using equations defined in Addison, Durrance and Shwietermann.

    Uses equation 7 in Addison, Durrance and Shwietermann to calculate the 
    total flux blocked per pixel solid angle, F_A, of a transiting 
    exoplanet. 
    '''


    ##########################################################################
    ############## List constants and determine basic parameters #############
    ##########################################################################

    # Define constants
    G = 6.673e-11 # Universal gravitational constant
    #m_ref = -1.5 # Apparent magnitude of Sirius
    #I_ref = (40*3.846e26)/(4*np.pi*(9.46e15)**2) # Intensity of Sirius. Units: W/m**2
    # N.B. Above is from I = L/4piR**2 where R = distance from Earth

    # Determine intensity of star (from apparent magnitude equation)
    #I_0 = I_ref*10**((2.5)*(m_ref-11.69))
    I_0 = 1


    ##########################################################################
    ########## Determine times over which transit has been observed ##########
    ##########################################################################

    # T is the time between observation start and the transit midpoint
    T = np.linspace(-(mid-start), +(end-mid), nobs)

    #########################################################################
    ################ Compute orbital velocities and angles ##################
    #########################################################################

    # Computes orbital velocity using two masses and the semi major-axis
    vorb = np.sqrt(G*(Mstar + Mplanet)/a)

    # Determine angular velocity where T = time from transit midpoint (s)
    #omega = vorb/a

    # Use above to calculate orbital phase angle
    #angle = T*omega
    #phase = np.linspace(-angle, angle, (nobs+1))


    ##########################################################################
    ############ Split exoplanet into arbitrary number of pixels #############
    ##########################################################################

    # Generates arbitrary coordinates in pixels and finds those that fit
    # within radius**2 (i.e. pythagorean). This is then used to calculate
    # these values as a fraction of the exoplanet radius.
    '''
    for x, y in product(np.linspace(-radius, radius, npix), repeat=2):
        if x ** 2 + y ** 2 <= radius**2:
            x_coord.append((x/radius)*Rplanet)
            y_coord.append((y/radius)*Rplanet)

    # Converts above lists to arrays for compatibility purposes
    x_coord = np.asarray(x_coord)
    y_coord = np.asarray(y_coord)
    '''

    x = np.linspace(-radius, +radius, (2*radius+1))
    y = np.linspace(-radius, +radius, (2*radius+1))

    coords = np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])

    x_coord = coords[:,0]
    y_coord = coords[:,1]

    x_exo, y_exo = [], []

    for i in range(len(coords)):
        if x_coord[i]**2 + y_coord[i]**2 <= radius**2:
            x_exo.append((x_coord[i]/radius)*Rplanet)
            y_exo.append((y_coord[i]/radius)*Rplanet)

    # Converts above lists to arrays for compatibility purposes
    x_exo = np.asarray(x_exo)
    y_exo = np.asarray(y_exo)

    ##########################################################################
    ################### Compute positions relative to star ###################
    ##########################################################################

    # Declare array to house values of x and y position
    X_pos = np.zeros(len(T))
    #Y_pos = np.zeros(len(phase))

    # Determine values of x position
    X_pos = vorb*T

    # Declare array to house flux values
    F_A = np.zeros(len(X_pos))

    # Determine y position (distance perpendicular to x position)
    #Y_pos = a*(np.arcsin(2*pi - radians(i)))    

    # Declare array to house distances from centre of star
    dpix = np.zeros((len(x_exo), len(X_pos)))

    ##########################################################################
    ############## Determine flux blocked per pixel solid angle ##############
    ##########################################################################
    
    # Determine position (postion from transit midpoint across stellar disk)
    for i in range(0, len(X_pos)):
        F_block = 0
        for j in range(0, len(x_exo)):
            
            #X_pos[i] = a*np.sin((phase[i]))
            #X_pos[i] = vorb*T[i]

            # Determine distance between centre of planetary disk and stellar disk
            # The numpy.sqrt allows sqrt of an array
            dpix_a = np.sqrt(((X_pos[i] + x_exo[j])**2) + y_exo[j]**2)
            #dpix[j] = dpix_a

            # Loop over intervals of d to calculate flux blocked. If the distance
            # to the exoplanet lies within the radius of the star, the flux blocked
            # is calculated using equation 7 of A, D and S paper. 
            # If it lies outside of the radius of the star, the flux blocked is assigned
            # a 0 value (as is in the array anyway). 
            if dpix_a <= Rstar:
                F_block += ((1-mu*(1-np.sqrt(1-(dpix_a/Rstar)**2)))/len(x_exo))

        F_A[i] = F_block

    ##########################################################################
    ##############  ##############
    ##########################################################################
       
    # Use simplified form of equation 7 to determine off-transit flux
    #F = np.ones(len(X_pos))*(I_0*(1-mu))
    F = np.ones(len(X_pos))

    # Subtract flux blocked from F to determine total flux per unit pixel solid angle 
    tot_F = F - F_A

    return tot_F, X_pos, F_A
