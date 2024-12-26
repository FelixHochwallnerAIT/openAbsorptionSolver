#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Felix Hochwallner, 2020
"""

#%% Import libraries

from CoolProp.CoolProp import PropsSI

#%% Functions

def hVap(T):
    """
    Returns latent heat of vaporization of water from temperature
    (using CoolProp)
    Input:
    T: Water temperature (degC)
    Output:
    hVap: Latent heat of vaporization of water (J/kg)
    """
    p       = PropsSI('P','T',max(T+273.15, 273.16),'Q',0.5,'Water')
    H_V     = PropsSI('H','P',p,'Q',1,'Water')
    H_L     = PropsSI('H','P',p,'Q',0,'Water')
    hVap  = H_V-H_L
    return hVap