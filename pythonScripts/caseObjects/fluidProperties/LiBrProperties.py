#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09.01.2020

@author: Felix Hochwallner
"""

#%% Import libraries

import numpy as np
from scipy import interpolate

#%% Functions

def rho(X,T):
    """
    Returns density of LiBr solution from concentration and temperature
    (Teja fit)
    Input:
    X: LiBr concentration (1)
    T: Solution temperature (degC)
    Output:
    rho: Density of solution (kg/m3)
    """
    A = np.array([  [ 1.097630E00,  0.071244E00,  2.21446E00],
                    [ 0.679620E-3, -1.482470E-3, -0.89696E-3],
                    [-0.035097E-6, -3.243120E-6,  4.97020E-6]])

    TT = (T+273.15)**np.array([0,1,2])
    XX = X**np.array([0,1,2])

    rho = np.sum(np.multiply(A.dot(XX),TT))*1e3
    return rho

def kappa(X,T):
    """
    Returns thermal conductivity of LiBr solution from concentration and 
    temperature (Patterson fit)
    Input:
    X: LiBr concentration (1)
    T: Solution temperature (degC)
    Output:
    kappa: Thermal conductivity of solution (W/mK)
    """
    A = np.array([  
        [   4.815196E-1, -2.217277E-3, -1.994141E-5,  
            3.727255E-7,  -2.489886E-9],
        [   1.858174E-3,  9.614755E-6, -1.139291E-6,  
            2.107608E-8,  -1.330532E-10],
        [   -7.923126E-6, -1.869392E-7,  1.408951E-8, 
            -2.740806E-10, 1.810818E-12]])

    TT = T**np.array([0,1,2])
    XX = (X*100)**np.array([0,1,2,3,4])

    kappa = np.sum(np.sum(np.multiply(A,np.outer(TT,XX))))*1.163
    return kappa

def cp(X,T):
    """
    Returns specific heat capacity of LiBr solution from concentration and 
    temperature (Feuerecker fit)
    Input:
    X: LiBr concentration (1)
    T: Solution temperature (degC)
    Output:
    cp: Specific heat capacity of solution (J/kgK)
    """
    a = np.array([-954.8, 47.7739, -1.59235, 2.09422E-2, -7.689E-5])
    b =  np.array([-3.293E-1, 4.076E-2, -1.36E-5,-7.1366E-6, 0])
    c =  np.array([7.4285E-3, -1.5144E-4, 1.3555E-6, 0, 0])
    d = -2.269E-6

    TT = T+273.15
    XX = (X*100)**np.array([0,1,2,3,4])

    cp = (np.sum(np.multiply(b,XX)) + 2*TT*np.sum(np.multiply(c,XX)) + 3*TT**2*d)*1e3
    return cp

def alpha(X,T):
    """
    Returns thermal diffusifity of LiBr solution from concentration and 
    temperature (calculated using kappa, rho and cp fits)
    Input:
    X: LiBr concentration (1)
    T: Solution temperature (degC)
    Output:
    alpha: Thermal diffusivity of solution (m2/s)
    """
    kappa_   = kappa(X,T)
    rho_     = rho(X,T)
    cp_      = cp(X,T)
    alpha    = kappa_/(rho_*cp_)
    return alpha

def mu(X,T):
    """
    Returns dynamic viscosity of LiBr solution from concentration and 
    temperature (Patterson fit)
    Input:
    X: LiBr concentration (1)
    T: Solution temperature (degC)
    Output:
    Mu: Dynamic viscosity of solution (Pas)
    """
    A = np.array([  
        [   1.488747E+0,  1.143975E-1, -1.278729E-2,  
            6.999985E-4, -1.638074E-5,  1.456348E-07],
        [   -4.164814E-2,  9.636832E-4, -5.981025E-5, 
            -1.282435E-7,  5.703002E-8, -9.842266E-10],
        [   3.404030E-4, -2.794515E-5,  2.580301E-6, 
            -9.737750E-8,  1.585609E-9, -7.922925E-12]])

    TT = T**np.array([0,1,2])
    XX = (X*100)**np.array([0,1,2,3,4,5])

    mu = np.sum(np.sum(np.multiply(A,np.outer(TT,XX))))/1e3
    return mu


def nu(X,T):
    """
    Returns kinematic viscosity of LiBr solution from concentration and 
    temperature (calculated using mu and rho fits)
    Input:
    X: LiBr concentration (1)
    T: Solution temperature (degC)
    Output:
    Nu: Kinematic viscosity of solution (m2/s)
    """
    nu = mu(X,T)/rho(X,T)
    return nu

def cwDiff(X,T):
    """
    Returns water diffusifity in LiBr solution from concentration and 
    temperature (Gierow fit)
    Input:
    X: LiBr concentration (1)
    T: Solution temperature (degC)
    Output:
    cwDiff: Water diffusivity in solution (m2/s)
    """
    T0 = 273.15

    Xv = np.array([ 1, 5, 11, 17, 23, 29.15, 35, 41, 47, 53, 59.27 ])/1e2
    Dv = np.array([ 1.321, 1.349, 1.44 , 1.539, 1.655, 1.739,
                    1.822, 1.826, 1.809, 1.488, 1.041 ])*1e-9

    eta25 = mu(X,25)

    # etaT = np.zeros(np.size(Xv))
    # DT = np.zeros(np.size(Xv))
    # for i in range(np.size(Xv)):
    #     etaT[i] = mu(Xv[i],T)
    #     DT[i]   = Dv[i]*eta25/(T0+25)*(T0+T)/etaT[i]
    
    f_D = interpolate.interp1d(
        Xv, Dv, 
        kind='cubic',
        bounds_error=False, 
        fill_value=(1.321e-9,1.041e-9))
    D25 = f_D(X)
    eta = mu(X,T)

    cwDiff =D25*eta25/(T0+25)*(T0+T)/eta
    return cwDiff