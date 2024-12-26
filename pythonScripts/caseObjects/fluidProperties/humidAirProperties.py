#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Felix Hochwallner, 2020
"""

#%% Import libraries

import numpy as np
from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI

#%% Functions

def Wcw(cw):
    """
    Returns humidity ratio from water concentration, i.e. specific humidity
    Input:
    cw: Water concentration, specific humidity (1) (kgW/kgA)
    Output:
    W: humidity ratio (1) (kgW/kgDA)
    """
    W = cw/(1-cw)
    return W

def cwW(W):
    """
    Returns water concentration, i.e. specific humidity from humidity ratio
    Input:
    W: humidity ratio (1) (kgW/kgDA)
    Output:
    cw: Water concentration, specific humidity (1) (kgW/kgA)
    """
    cw = W/(1+W)
    return cw

def cwPw(p,pw):
    """
    Returns water concentration, i.e. specific humidity, from system and 
    water vapor pressure
    Input:
    p: Air pressure (Pa)
    pw: Water vapor pressure (Pa)
    Output:
    cw: Water concentration, specific humidity (1) (kgW/kgA)
    """
    beta = 0.6225
    cw = beta*pw/(p-(1-beta)*pw)
    return cw

def cwRH(p,T,RH):
    """
    Returns water concentration, i.e. specific humidity, from humid air, 
    given the air pressure, temperature and relative humidity.
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Air relative humidity (1)
    Output:
    cw: Water concentration, specific humidity (1) (kgW/kgA)
    """
    W = HAPropsSI('W', 'P', p, 'T', T+273.15, 'RH', RH)
    cw = cwW(W)
    return cw

def RHpTcw(p,T,cw):
    """
    Returns relative humidity of humid air, given the air
    pressure, temperature and specific humidity.
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    cw: Specific humidity (1) (kgW/kgA)
    Output:
    RH: Relative humidity (1)
    """
    W = Wcw(cw)
    RH = HAPropsSI('RH', 'P', p, 'T', T+273.15, 'W', W)
    return RH

def pwRH(p,T,RH):
    """
    Returns water vapor pressure from humid air, given the air
    pressure, temperature and relative humidity.
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Air relative humidity (1)
    Output:
    pwRH: Water vapor pressure (Pa)
    """
    pw = HAPropsSI('P_w', 'P', p, 'T', T+273.15, 'RH', RH)
    return pw
    
def pwC(p,T,C):
    """
    Returns water vapor pressure from humid air, given the air
    pressure, temperature and specific humidity.
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    C: Air specific humidity (1)
    Output:
    pw: Water vapor pressure (Pa)
    """
    pw = HAPropsSI('P_w', 'P', p, 'T', T+273.15, 'HumRat', C)
    return pw

def rho(p, T, RH):
    """
    Returns density of air from pressure, temperature and relative humidity 
    (using CoolProp)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Relative humidity (%)
    Output:
    rho: Air density (kg/m3)
    """
    rho = 1/HAPropsSI('Vha', 'P', p, 'T', T+273.15, 'RH', RH)
    return rho

def kappa(p, T, RH):
    """
    Returns thermal conductivity of air from pressure, temperature and relative
    humidity (using CoolProp)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Relative humidity (%)
    Output:
    kappa: Air thermal conductivity (W/mK)
    """
    kappa = HAPropsSI('K', 'P', p, 'T', T+273.15, 'RH', RH)
    return kappa

def cp(p, T, RH):
    """
    Returns specific heat capacity of air from pressure, temperature and 
    relative humidity (using CoolProp)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Relative humidity (%)
    Output:
    cp: Specific heat capacity of air (J/kgK)
    """
    cp = HAPropsSI('cp_ha', 'P', p, 'T', T+273.15, 'RH', RH)
    return cp

def cwDiff(p, T, which='VDI'):
    """
    Returns water diffusifity in air from pressure and temperature 
    (constant)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    Output:
    cwDiff: Water diffusivity in air (m2/s)
    """
    if which=='Yonko':
        cwDiff = 2.302*(p/0.98e5)*((T+273.15)/256)**1.81*10**(-5)
            # Yonko (see Schmitz paper)
    elif which=='Schirmer':
        cwDiff = 0.083/3600 * ((T+273.15)/273)**(1.81) 
            # Schirmer (1938)
    elif which=='VDI':
        cwDiff = 2.252/p * ((T+273.15)/273)**(1.81) 
            # VDI (2013)
    else:
        cwDiff = 2.5e-5
    return cwDiff

def alpha(p, T, RH):
    """
    Returns thermal diffusifity of air from pressure, temperature and 
    relative humidity (using CoolProp)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Relative humidity (%)
    Output:
    alpha: Thermal diffusivity of air (m2/s)
    """
    kappa_   = kappa(p, T, RH)
    rho_     = rho(p, T, RH)
    cp_      = cp(p, T, RH)
    alpha    = kappa_/(rho_*cp_)
    return alpha

def mu(p, T, RH):
    """
    Returns dynamic viscosity of air from pressure, temperature and 
    relative humidity (using CoolProp)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Relative humidity (%)
    Output:
    Nu: Dynamic viscosity of air (Pas)
    """
    mu = HAPropsSI('mu', 'P', p, 'T', T+273.15, 'RH', RH)
    return mu

def nu(p, T, RH):
    """
    Returns kinematic viscosity of air from pressure, temperature and 
    relative humidity (using CoolProp)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Relative humidity (%)
    Output:
    Nu: Kinematic viscosity of air (m2/s)
    """
    nu = mu(p, T, RH)/rho(p, T, RH)
    return nu

def Pr(p, T, RH):
    """
    Returns Prandtl number of air from pressure, temperature and 
    relative humidity (using CoolProp)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Relative humidity (%)
    Output:
    Pr: Prandtl number of air (1)
    """
    Pr = nu(p, T, RH)/alpha(p, T, RH)
    return Pr

def beta():
    """
    Returns the molar mass ratio of water to air from pressure, temperature 
    and relative humidity (using CoolProp)
    Input:
    p: Air pressure (Pa)
    T: Air temperature (degC)
    RH: Relative humidity (%)
    Output:
    beta: Prandtl number of air (1)
    """
    M_water = PropsSI('M','T',273,'P',1e5,'Water')
    M_air = PropsSI('M','T',273,'P',1e5,'Air')
    beta = M_water / M_air
    return beta
    