#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:17:21 2020

@author: Felix Hochwallner

Class definition initializationSettings, for multiRegionFallingFilm OpenFoam
class cases.
"""

class airSolutionInitializationSettings:
    """
    Stores the initialization settings of an OpenFoam case.
    """
    def __init__(self, air_init, solution_init):
        """
        Constructor class
        Input:
        air_init: airInitializationSettings object
        solution_init: solutionInitializationSettings object
        """
        self.air = air_init
        self.solution = solution_init


class airInitializationSettings:
    """
    Stores the air initialization settings.
    """
    def __init__(self, v_in=2, T_in=20, C_in=14e-3):
        """
        Constructor class
        Input:
        v_in: Air velocity in x-direction at inlet (m/s)
        T_in: Air temperature at inlet (degC)
        C_in: Water concentration at inlet (1)
        """
        self.v_in = v_in
        self.T_in = T_in
        self.C_in = C_in

class solutionInitializationSettings:
    """
    Stores the air initialization settings.
    """
    def __init__(self, T_in=20, C_in=0.2, T_wall=20, mDot=1.768e-3):
        """
        Constructor class
        Input:
        T_in: Solution temperature at inlet (degC)
        C_in: Water concentration at inlet (1)
        T_wall: Wall temperature (degC)
        """
        self.T_in = T_in
        self.C_in = C_in
        self.T_wall = T_wall
        self.mDot = mDot
        