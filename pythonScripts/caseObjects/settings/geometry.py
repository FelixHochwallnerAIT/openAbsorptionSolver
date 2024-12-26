#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:17:21 2020

@author: Felix Hochwallner

Class definition geometrySettings, for OpenFoam case classes.
"""

import numpy as np
from PyFoam.Basics.DataStructures import Vector

class airWireGeometrySettings:
    """
    Stores the geometry settings of a wire air flow OpenFoam case.
    """
    def __init__(self, 
        l=0.1, b=2e-3, h=0.7, d_wall=5e-4, d_wire=2e-3, n_wires=1,
        r=[0.25], k1=5.0, k2=20.0,
        x_inlet=10e-3, x_outlet=10e-3, 
        grading_inlet_outlet=False, grading_inlet_outlet_ratio=1.2,
        ny_air=40, ny_sol=10, nz=40):
        """
        Constructor class
        Input:
        l: Length of air channel (m)
        b: Width of air channel (m)
        h: Height of air channel (m)
        etc.
        """
        self.type = 'wireWithInletOutlet'
        # Geometric dimensions (m)
        self.l = l
        self.b = b
        self.h = h
        self.d_wall = d_wall
        self.l_total = self.l

        # Wire data
        self.d_wire = d_wire    # Wire diameter
        self.n_wires = n_wires  # Number of wires
        self.k1 = k1            # Entrance length
        self.k2 = k2            # Exit length

        # Inlet and outlet
        self.x_inlet = x_inlet
        self.x_outlet = x_outlet
        self.grading_inlet_outlet = grading_inlet_outlet
        self.grading_inlet_outlet_ratio = grading_inlet_outlet_ratio

        # Set amount of cells
        self.ny_air = ny_air
        self.ny_sol = ny_sol
        self.nz = nz
        # dx = (self.b-2*d_Nu-d_wall) / self.ny_air
        # self.dx = dx
        self.ny1 = int(self.ny_air / 4)
        assert ny_air % 4 == 0, "ERROR: ny_air must be a multiple of 4."

        self.r = r              # Ratio where wire are located
        assert len(self.r) == self.n_wires, (
            "ERROR: Number of wires and ratios do not match.")

class solutionWireGeometrySettings:
    """
    Stores the geometry settings of a wire solution flow OpenFoam case.
    """
    def __init__(self, 
        l=0.1, b=2e-3, h=0.7, d_wall=5e-4, d_wire=2e-3, n_wires=1,
        r=[0.25], k1=5.0, k2=20.0,
        x_inlet=10e-3, x_outlet=10e-3,
        grading_inlet_outlet=False, grading_inlet_outlet_ratio=1.2,
        ny_air=40, ny_sol=10, nz=40,
        grading=1):
        """
        Constructor class
        Input:
        l: Length of air channel (m)
        b: Width of air channel (m)
        h: Height of air channel (m)
        etc.
        """
        self.type = 'wireWithInletOutlet'
        # Geometric dimensions (m)
        self.l = l
        self.b = b
        self.h = h
        self.d_wall = d_wall
        self.l_total = self.l

        # Wire data
        self.d_wire = d_wire    # Wire diameter
        self.n_wires = n_wires  # Number of wires
        self.k1 = k1            # Entrance length
        self.k2 = k2            # Exit length

        # Inlet and outlet
        self.x_inlet = x_inlet
        self.x_outlet = x_outlet
        self.grading_inlet_outlet = grading_inlet_outlet
        self.grading_inlet_outlet_ratio = grading_inlet_outlet_ratio

        # Set amount of cells
        self.ny_air = ny_air
        self.ny_sol = ny_sol
        self.nz = nz
        self.grading = grading
        # dx = self.b / self.ny_air
        self.ny1 = int(self.ny_air / 4)

        self.r = r              # Ratio where wire are located
        assert len(self.r) == self.n_wires, (
            'ERROR: Number of wires and ratios do not match.')

class wireGeometrySettings:
    """
    Stores the geometry settings of a wire OpenFoam case.
    """
    def __init__(self, airGeometry, solutionGeometry):
        """
        Constructor class
        Input:
        airGeometry: airGeometrySettings object
        solutionGeometry: solutionGeometrySettings object
        """
        self.air = airGeometry
        self.solution = solutionGeometry
        if not (airGeometry.type == solutionGeometry.type):
            print("ERROR: Different mesh types for air and solution:")
            print("air: %s; solution: %s" %(self.air.type, self.solution.type))
            raise
        self.type = solutionGeometry.type
        # Set properties to global attribute
        self.l = solutionGeometry.l
        self.b = solutionGeometry.b
        self.h = solutionGeometry.h
        self.d_wall = solutionGeometry.d_wall
        self.l_total = self.l
        # Wire data
        self.d_wire = solutionGeometry.d_wire    # Wire diameter
        self.n_wires = solutionGeometry.n_wires  # Number of wires
        self.k1 = solutionGeometry.k1            # Entrance length
        self.k2 = solutionGeometry.k2            # Exit length
        self.r = solutionGeometry.r
        # Inlet and outlet
        self.x_inlet = solutionGeometry.x_inlet
        self.x_outlet = solutionGeometry.x_outlet
        self.grading_inlet_outlet = solutionGeometry.grading_inlet_outlet
        self.grading_inlet_outlet_ratio = (
            solutionGeometry.grading_inlet_outlet_ratio)
        # Set amount of cells
        self.ny_air = solutionGeometry.ny_air
        self.ny_sol = solutionGeometry.ny_sol
        self.nz = solutionGeometry.nz
        # dx = (self.b-2*d_Nu-d_wall) / self.ny_air
        # self.dx = dx
        self.ny1 = int(self.ny_air / 4)

class airStackedArrangementGeometrySettings:
    """
    Stores the geometry settings of stacked arrangement air 
    flow OpenFoam case.
    """
    def __init__(self, l=0.1, b=4e-3, h=0.7, d_wall=5e-4, n_steps=1,
        ny_air=30, ny_sol=10, nz=50, x_gap=10e-3, 
        x_inlet=10e-3, x_outlet=10e-3, 
        grading_inlet_outlet=False, grading_inlet_outlet_ratio=1.2):
        """
        Constructor class
        Input:
        l: Length of air channel (m)
        b: Width of air channel (m)
        h: Height of air channel (m)
        d_wall: Width of wall (m)
        lr_inlet: Ratio of inlet zone (1)
        rnx_inlet: Ratio of amount of cells in inlet zone (1)
        nx: Amount of cells in x-dir (length) of air channel (1)
        ny: Amount of cells in y-dir (width) of air channel (1)
        nz: Amount of cells in z-dir (height) of air channel (1)
        """
        self.type = 'stackedArrangementWithGapAndInletOutlet'
        # Geometric dimensions
        self.l = l
        self.b = b
        self.h = h
        self.d_wall = d_wall
        self.x_gap = x_gap
        self.x_inlet = x_inlet
        self.x_outlet = x_outlet
        self.grading_inlet_outlet = grading_inlet_outlet
        self.grading_inlet_outlet_ratio = grading_inlet_outlet_ratio
        self.n_steps = n_steps
        self.l_total = (
            self.l * self.n_steps
            + self.x_gap * (self.n_steps-1))
        # Amount of cells
        self.ny_air = ny_air
        self.ny_sol = ny_sol
        self.nz = nz

class solutionStackedArrangementGeometrySettings:
    """
    Stores the geometry settings of stacked arrangement solution 
    flow OpenFoam case.
    """
    def __init__(self, l=0.1, b=4e-3, h=0.7, d_wall=5e-4, n_steps=1,
        ny_air=30, ny_sol=10, nz=50, grading=1, x_gap=10e-3,
        x_inlet=10e-3, x_outlet=10e-3, 
        grading_inlet_outlet=False, grading_inlet_outlet_ratio=1.2):
        """
        Constructor class
        Input:
        l: Length of air channel (m)
        b: Width of air channel (m)
        h: Height of air channel (m)
        d_wall: Width of wall (m)
        lr_inlet: Ratio of inlet zone (1)
        rnx_inlet: Ratio of amount of cells in inlet zone (1)
        nx: Amount of cells in x-dir (length) of air channel (1)
        ny: Amount of cells in y-dir (width) of air channel (1)
        nz: Amount of cells in z-dir (height) of air channel (1)
        """
        self.type = 'stackedArrangementWithGapAndInletOutlet'
        # Geometric dimensions
        self.l = l
        self.b = b
        self.h = h
        self.d_wall = d_wall
        self.x_gap = x_gap
        self.x_inlet = x_inlet
        self.x_outlet = x_outlet
        self.grading_inlet_outlet = grading_inlet_outlet
        self.grading_inlet_outlet_ratio = grading_inlet_outlet_ratio
        self.n_steps = n_steps
        self.l_total = self.l
        # Amount of cells
        self.ny_air = ny_air
        self.ny_sol = ny_sol
        self.nz = nz
        self.grading = grading

class stackedArrangementGeometrySettings:
    """
    Stores the geometry settings of a stacked arrangement OpenFoam case.
    """
    def __init__(self, airGeometry, solutionGeometry, 
        overwrite_check=False, calc_ny_sol=False):
        """
        Constructor class
        Input:
        airGeometry: airGeometrySettings object
        solutionGeometry: solutionGeometrySettings object
        """
        self.overwrite_check = overwrite_check
        self.calc_ny_sol = calc_ny_sol
        self.air = airGeometry
        self.solution = solutionGeometry
        if not (airGeometry.type == solutionGeometry.type):
            print("ERROR: Different mesh types for air and solution:")
            print("air: %s; solution: %s" %(self.air.type, self.solution.type))
            raise
        self.type = solutionGeometry.type
        # Set properties to global attribute
        self.l = solutionGeometry.l
        self.b = solutionGeometry.b
        self.h = solutionGeometry.h
        self.d_wall = solutionGeometry.d_wall 
        self.x_gap = solutionGeometry.x_gap
        self.x_inlet = solutionGeometry.x_inlet
        self.x_outlet = solutionGeometry.x_outlet
        self.grading_inlet_outlet = solutionGeometry.grading_inlet_outlet
        self.grading_inlet_outlet_ratio = (
            solutionGeometry.grading_inlet_outlet_ratio)
        self.n_steps = solutionGeometry.n_steps
        self.l_total = airGeometry.l_total
        # Amount of cells
        self.ny_air = solutionGeometry.ny_air
        self.ny_sol = solutionGeometry.ny_sol
        self.nz = solutionGeometry.nz
