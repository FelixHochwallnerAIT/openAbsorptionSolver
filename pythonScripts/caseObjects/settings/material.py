#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:17:21 2020

@author: Felix Hochwallner

Class definitions for material settings, for OpenFoam case classes.
"""

class airSolutionMaterialSettings:
    """
    Stores the air and solution material settings of a multi region 
    OpenFoam case.
    """
    def __init__(self, airMaterial, solutionMaterial):
        """
        Constructor class
        Input:
        airMaterial: airMaterialSettings object
        solutionMaterial: solutionMaterialSettings object
        """
        self.air = airMaterial
        self.solution = solutionMaterial

class airMaterialSettings:
    """
    Stores the air material settings.
    """
    def __init__(self, rho=1.2, cp=1e3, mu=18.6e-6, Pr=0.707, kappa=26.24e-3, 
            cwDiff=2.5e-5, p0=1e5, beta=0.6225, deltaH=2.44e6):
        """
        Constructor class
        """
        self.rho = rho # kg/m3
        self.cp = cp # J/kgK
        self.mu = mu # Pas
        self.Pr =Pr # 1
        self.kappa = kappa # W/mK  kappa = lambda
        self.cwDiff = cwDiff # m2/s
        self.p0 = p0 # Pa
        self.beta = beta # 1
        self.deltaH = deltaH # J/kg
        self.nu = self.mu/self.rho # m2/s
        self.alpha = self.kappa/(self.rho*self.cp) # m2/s

class solutionFitSettings:
    """
    Stores the solution fit material settings.
    """
    def __init__(self,medium='default'):
        """
        Constructor class
        """
        if (medium == 'LiBr' or medium == 'default'):
            # Derived from Patterson (1984)
            self.k0A = 1.33790476e+01
            self.k1A = 1.27189268e+01
            self.k2A = -1.21799427e+01
            self.k3A = 0
            self.k0B = 3.95458129e+03
            self.k1B = 8.86537780e+03
            self.k2B = -1.90328499e+03
            self.k3B = 0
            self.k0C = 1.02517469e+02
            self.k1C = 3.50866545e+02
            self.k2C = 0
            self.k3C = 0
            self.D   = 0
            if medium=='default':
                print("LiBr as default desiccant selected.")
        # ADD NEW MEDIUM HERE
        else:
            print("ERROR: No known desiccant %s selected." %medium)
            raise

class solutionMaterialSettings:
    """
    Stores the solution material settings.
    """
    def __init__(self, medium='default'):
        """
        Constructor class
        Input:
        medium: Name of medium (string)
        """
        
        if (medium == 'LiBr' or medium == 'default'):
            self.medium = 'LiBr'
            self.rho = 1.5648e3                         # kg/m3
            self.cp = 2.1127e3                          # J/kgK
            self.mu = 0.0035                            # Pas
            self.kappa = 0.4432                         # W/mK  kappa = lambda
            self.cwDiff = 1.6841e-9                     # m2/s
            self.p0 = 1e5                               # Pa
            self.hs = 2.44e6                            # J/kg
            self.nu = self.mu/self.rho                  # m2/s
            self.alpha = self.kappa/(self.rho*self.cp)  # m2/s
            self.Pr = self.mu/self.rho/self.alpha       # 1
        # ADD NEW MEDIUM HERE
        else:
            print("ERROR: No known desiccant %s selected." %medium)
            raise
        
        # Initalize fit but deactivate it
        self.fit = solutionFitSettings(medium)
        self.flagFit = 0

    def setFit(self):
        """
        Setting solution fit and activating fit
        Input:
        fit: Solution fit (solutionFitSettings object)
        """
        self.flagFit = 1
