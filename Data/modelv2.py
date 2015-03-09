'''
Project: Exoplanetary Detection and Characterisation
Supervisor: Dr. Edward L Gomez

Author: Tomas James

Script: Prototype for modelling exoplanetary lightcurves
'''

import numpy as np
from numpy import asarray
from matplotlib.pyplot import *
from itertools import product

def model(Mstar, Rstar, Mplanet, Rplanet, radius, obs_end, trans_mid, nobs, a, i, mu, app_mag):

    '''
    Takes a planetary system consisting of an exoplanet of mass M_planet
    orbiting a star of mass M_star and radius R_star and models a light-
    curve using equations defined in Addison, Durrance and Shwietermann.

    Uses equation 7 in Addison, Durrance and Shwietermann to calculate the 
    total flux blocked per pixel solid angle, F_A, of a transiting 
    exoplanet. To do this, the duration of the transit (dur) along with the
    semi-major axis of the planet-star system  (a)and orbital inclination (i) are
    required. Furthermore the limb-darkening coefficient (mu) is needed, along
    with the apparent magnitude of the star under consideration (app_mag). 
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
    #I_0 = I_ref*10**((2.5)*(m_ref-app_mag))
    I_0 = 1


    ##########################################################################
    ########## Determine times over which transit has been observed ##########
    ##########################################################################

    # T is the time between observation start and the transit midpoint
    T = (obs_end - trans_mid)


#   #########################################################################
    ################# Compute orbital velocities and angles #################

    # Computes orbital velocity using two masses and the semi major-axis
    vorb = np.sqrt(G*(Mstar + Mplanet)/a)

    # Determine angular velocity where T = time from transit midpoint (s)
    omega = vorb/a

    # Use above to calculate orbital phase angle
    angle = T*omega
    phase = np.linspace(-angle, angle, (nobs+1))


    ##########################################################################
    ############ Split exoplanet into arbitrary number of pixels #############
    ##########################################################################

    # Declare empty list to store values of coordinates
    x_coord, y_coord = [], []

    # Generates arbitrary coordinates in pixels and finds those that fit
    # within radius**2 (i.e. pythagorean). This is then used to calculate
    # these values as a fraction of the exoplanet radius.
    for x, y in product(range(-radius, radius), repeat=2):
        if x ** 2 + y ** 2 <= radius**2:
            x_coord.append(x*(Rplanet/radius))
            y_coord.append(y*(Rplanet/radius))

    # Converts above lists to arrays for compatibility purposes
    x_coord = asarray(x_coord)
    y_coord = asarray(y_coord)


    ##########################################################################
    ################### Compute positions relative to star ###################
    ##########################################################################

    # Declare array to house values of x and y position
    X_pos = np.zeros(len(phase))
    Y_pos = np.zeros(len(phase))

    # Declare array to house flux values
    F_A = np.zeros(len(X_pos))

    # Determine y position (distance perpendicular to x position)
    #Y_pos = a*(np.arcsin(2*pi - radians(i)))    

    # Declare array to house distances from centre of star
    D = np.zeros(len(phase))

    # Determine position (postion from transit midpoint across stellar disk)
    for i in range(0, len(phase)):
        for j in range(0, len(x_coord)):
            
            X_pos[i] = a*np.sin((phase[i]))

            # Determine distance between centre of planetary disk and stellar disk
            # The numpy.sqrt allows sqrt of an array
            D[i] = np.sqrt(((X_pos[i] + x_coord[j])**2) + y_coord[j]**2)

            # Loop over intervals of d to calculate flux blocked. If the distance
            # to the exoplanet lies within the radius of the star, the flux blocked
            # is calculated using equation 7 of A, D and S paper. 
            # If it lies outside of the radius of the star, the flux blocked is assigned
            # a 0 value. 
            if D[i] <= Rstar:
                F_A[i] += I_0*(1-mu*(1-np.sqrt(1-(D[i]/Rstar)**2)))
            else:
                F_A[i] = 0
    

    ##########################################################################
    ############## Determine flux blocked per pixel solid angle ##############
    ##########################################################################
       
    # Use simplified form of equation 7 to determine off-transit flux
    #F = np.ones(len(X_pos))*(I_0*(1-mu))
    F = np.ones(len(X_pos))

    # Subtract flux blocked from F to determine total flux per unit pixel solid angle 
    tot_F = F - F_A

    return tot_F, X_pos, D
