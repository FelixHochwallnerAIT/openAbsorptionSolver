#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
15.05.2023

Felix Hochwallner

Class definition wire simulations, to view, edit and simulate 
prescribed wire cases.
"""

# Import standard python modules
import os
import numpy as np
import sys
import copy
from scipy.optimize import newton

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.ParsedBlockMeshDict import ParsedBlockMeshDict
from PyFoam.Basics.DataStructures import Vector as OFVector
from PyFoam.Basics.DataStructures import Field as OFField

from Ofpp import parse_internal_field, parse_field_all, parse_boundary_field

from .fluidProperties import humidAirProperties
from .baseCaseOpenFoam import baseCaseOpenFoam

class wire(baseCaseOpenFoam):
    """
    Stores an OpenFoam wire case.
    """
    def __init__(self, solver, geometry, material, initialization,
        case=''):
        """
        Extend constructor class
        """
        baseCaseOpenFoam.__init__(
            self, 
            solver, 
            geometry, 
            material, 
            initialization)
        self.classtype = 'wire'
        if not self.solver.multiRegion:
            print("ERROR: Wire case only supports multiRegion cases.")
            raise
        self.solver.setTransient()
        self.setFilmThickness()
        self.initalize(case=case)

    def calculateVertices(self):
        """
        Calculates the vertices of the blockMesh
        """
        geometry = self.geometry
        # Reduce falling film thickness by wall and two falling films
        geometry.b = geometry.b-2*geometry.d_Nu-geometry.d_wall
        assert geometry.b > 0, ("Air gap is too small to fit wall "
            "and two falling films")
        # Calculate dx 
        dx = (geometry.b) / geometry.ny_air
        geometry.dx = dx
        # Calculate auxiliary lengths
        geometry.y0 = geometry.d_wall/2.0
        geometry.y1 = (geometry.b/2-geometry.d_wire/2.0)/2.0
        geometry.x1 = ((2.0)**(1/2)
            *(geometry.b+geometry.d_wire)/8.0)
        geometry.y2 = ((
            4.0*geometry.b
            -(2.0)**(1/2)*geometry.b
            -(2.0)**(1/2)*geometry.d_wire)/8.0)   
            
        # Lengths in x-dir
        geometry.x = [
            0]
        for i in range(geometry.n_wires):
            geometry.x.append(
                geometry.r[i]*geometry.l
                -geometry.x1
                -geometry.y2*geometry.k1)
            geometry.x.append(
                geometry.r[i]*geometry.l
                -geometry.d_wire/2.0
                -geometry.y1)
            geometry.x.append(
                geometry.r[i]*geometry.l
                -geometry.x1)
            geometry.x.append(
                geometry.r[i]*geometry.l
                -geometry.d_wire/2.0)
            geometry.x.append(
                geometry.r[i]*geometry.l
                -geometry.d_wire/2.0*(2.0)**(1/2)/2.0)
            geometry.x.append(
                geometry.r[i]*geometry.l)
            geometry.x.append(
                geometry.r[i]*geometry.l
                +geometry.d_wire/2.0*(2.0)**(1/2)/2.0)
            geometry.x.append(
                geometry.r[i]*geometry.l
                +geometry.d_wire/2.0)
            geometry.x.append(
                geometry.r[i]*geometry.l
                +geometry.x1)
            geometry.x.append(
                geometry.r[i]*geometry.l
                +geometry.d_wire/2.0+geometry.y1)
            geometry.x.append(
                geometry.r[i]*geometry.l
                +geometry.x1
                +geometry.y2*geometry.k2)
        geometry.x.append(geometry.l)

        # Lengths in y-dir
        geometry.y = [
            geometry.y0-geometry.d_Nu,
            geometry.y0,
            geometry.y0+geometry.y1,
            geometry.y0+geometry.y2,
            geometry.y0+2*geometry.y1,
            geometry.y0+geometry.b/2-(2.0)**(1/2)*geometry.d_wire/4.0,
            geometry.y0+geometry.b/2,
            geometry.y0+geometry.b/2+(2.0)**(1/2)*geometry.d_wire/4.0,
            geometry.y0+geometry.b/2+geometry.d_wire/2.0,
            geometry.y0+geometry.b/2+geometry.b/2-geometry.y2,
            geometry.y0+geometry.b/2+geometry.d_wire/2.0+geometry.y1,
            geometry.y0+geometry.b,
            geometry.y0+geometry.b+geometry.d_Nu]
        geometry.yh = [
            geometry.y0+geometry.b/4.0,
            geometry.y0+geometry.b*3.0/4.0]
        # Lengths in z-dir
        geometry.z = [
            0,
            geometry.h]

        # Rounding points
        # x-dir
        geometry.xp = []
        for i in range(geometry.n_wires):
            geometry.xp.append(
                geometry.r[i]*geometry.l
                -(3.0)**(1/2)*(geometry.d_wire+geometry.b)/8.0)
            geometry.xp.append(
                geometry.r[i]*geometry.l
                -(3.0)**(1/2)*geometry.d_wire/4.0)
            geometry.xp.append(
                geometry.r[i]*geometry.l
                -(geometry.d_wire+geometry.b)/8.0)
            geometry.xp.append(
                geometry.r[i]*geometry.l
                -geometry.d_wire/4.0)
            geometry.xp.append(
                geometry.r[i]*geometry.l
                +geometry.d_wire/4.0)
            geometry.xp.append(
                geometry.r[i]*geometry.l
                +(geometry.d_wire+geometry.b)/8.0)
            geometry.xp.append(
                geometry.r[i]*geometry.l
                +(3.0)**(1/2)*geometry.d_wire/4.0)
            geometry.xp.append(
                geometry.r[i]*geometry.l
                +(3.0)**(1/2)*(geometry.d_wire+geometry.b)/8.0)
        # y-dir
        geometry.yp = [
            geometry.y0+geometry.b/2-(3.0)**(1/2)
                *(geometry.d_wire+geometry.b)/8.0,
            geometry.y0+geometry.b/2-(3.0)**(1/2)*geometry.d_wire/4.0,
            geometry.y0+geometry.b/2-(geometry.d_wire+geometry.b)/8.0,
            geometry.y0+geometry.b/2-geometry.d_wire/4.0,
            geometry.y0+geometry.b/2+geometry.d_wire/4.0,
            geometry.y0+geometry.b/2+(geometry.d_wire+geometry.b)/8.0,
            geometry.y0+geometry.b/2+(3.0)**(1/2)*geometry.d_wire/4.0,
            geometry.y0+geometry.b/2+(3.0)**(1/2)
                *(geometry.d_wire+geometry.b)/8.0]
        geometry.nx = [
            int((geometry.x[1] - geometry.x[0]) / dx),
            int(geometry.y2*geometry.k1 / dx),
            int(geometry.ny1),
            int(geometry.ny1),
            int(geometry.y2*geometry.k2 / dx)]

        def find_grading(x, delta_s, ratio_max=1.2):
            # Find n
            n = 0
            x_temp = x
            delta_s_temp = delta_s
            while x_temp > 0:
                x_temp -= delta_s_temp
                delta_s_temp *= ratio_max
                n += 1
            # Find ratio
            def ratio_func(ratio):
                ret = 0
                for i in range(n):
                    ret += delta_s * (ratio**i)
                ret -= x
                return ret
            try:
                ratio = newton(ratio_func, 1.06, tol=1e-8, maxiter=100000)
            except:
                ratio = 1
                n = int(round(x/delta_s))
            # ret = ratio_func(ratio)
            # Find grading
            grading = ratio**(n-1)
            return n, grading
        if geometry.grading_inlet_outlet:
            geometry.nx_inlet, geometry.grading_inlet = find_grading(
                geometry.x_inlet, 
                dx,
                geometry.grading_inlet_outlet_ratio)
            geometry.grading_inlet = 1/geometry.grading_inlet
            geometry.nx_outlet, geometry.grading_outlet = find_grading(
                geometry.x_outlet, 
                dx,
                geometry.grading_inlet_outlet_ratio)
        else:
            geometry.nx_inlet = int(geometry.x_inlet / dx)
            geometry.nx_outlet = int(geometry.x_outlet / dx)
        geometry.ny_wall = int(np.ceil(geometry.d_wall/2 / dx))
        geometry.ny_sol_in_out = int(np.ceil(geometry.d_Nu / dx))
        for i in range(geometry.n_wires-1):
            geometry.nx.append(int((
                geometry.r[i+1]*geometry.l 
                - geometry.x1 
                - geometry.y2*geometry.k1 
                - geometry.r[i]*geometry.l 
                - geometry.x1 
                - geometry.y2*geometry.k2) 
                / dx))
        if geometry.n_wires > 0:
            geometry.nx.append(int((
                geometry.l 
                - geometry.r[-1]*geometry.l 
                - geometry.x1 
                - geometry.y2*geometry.k2)
                / dx))
        # Number of cells
        geometry.n_cells = {'total': 0}
        for region in self.solver.regions:
            if region == 'air':
                # Blocks 0 to 3
                geometry.n_cells[region] = int(
                    4 * geometry.nx[0] * geometry.ny1 * geometry.nz)
                # Remaining blocks
                for i in range(geometry.n_wires):
                    geometry.n_cells[region] += int(
                        4 * geometry.nx[1] * geometry.ny1 * geometry.nz)
                    geometry.n_cells[region] += int(
                        6 * geometry.nx[2] * geometry.ny1 * geometry.nz)
                    geometry.n_cells[region] += int(
                        6 * geometry.nx[3] * geometry.ny1 * geometry.nz)
                    geometry.n_cells[region] += int(
                        4 * geometry.nx[4] * geometry.ny1 * geometry.nz)
                    geometry.n_cells[region] += int(
                        4 * geometry.nx[5+i] * geometry.ny1 * geometry.nz)
                # outlet
                geometry.n_cells[region] += int(
                    4 * geometry.nx_outlet * geometry.ny1 * geometry.nz)
                geometry.n_cells[region] += int(
                    2 * geometry.nx_outlet 
                    * geometry.ny_sol_in_out 
                    * geometry.nz)
                geometry.n_cells[region] += int(
                    2 * geometry.nx_outlet * geometry.ny_wall * geometry.nz)
                # inlet
                geometry.n_cells[region] += int(
                    4 * geometry.nx_inlet * geometry.ny1 * geometry.nz)
                geometry.n_cells[region] += int(
                    2 * geometry.nx_inlet 
                    * geometry.ny_sol_in_out 
                    * geometry.nz)
                geometry.n_cells[region] += int(
                    2 * geometry.nx_inlet * geometry.ny_wall * geometry.nz)
            elif region == 'solution1' or region == 'solution2':
                # Blocks 0 to 3
                geometry.n_cells[region] = int(
                    geometry.nx[0] * geometry.ny_sol * geometry.nz)
                # Remaining blocks
                for i in range(geometry.n_wires):
                    geometry.n_cells[region] += int(
                        geometry.nx[1] * geometry.ny_sol * geometry.nz)
                    geometry.n_cells[region] += int(
                        geometry.nx[2] * geometry.ny_sol * geometry.nz)
                    geometry.n_cells[region] += int(
                        geometry.nx[3] * geometry.ny_sol * geometry.nz)
                    geometry.n_cells[region] += int(
                        geometry.nx[4] * geometry.ny_sol * geometry.nz)
                    geometry.n_cells[region] += int(
                        geometry.nx[5+i] * geometry.ny_sol * geometry.nz)
            else:
                print(f'ERROR: Region {region} not defined.')
                sys.exit(1)
            geometry.n_cells['total'] += geometry.n_cells[region]

        # Create vertices
        vertices = [
            OFVector(geometry.x[0], geometry.y[0], geometry.z[0]),
            OFVector(geometry.x[0], geometry.y[0], geometry.z[1]),
            OFVector(geometry.x[0], geometry.y[1], geometry.z[0]),
            OFVector(geometry.x[0], geometry.y[1], geometry.z[1]),
            OFVector(geometry.x[0], geometry.yh[0], geometry.z[0]),
            OFVector(geometry.x[0], geometry.yh[0], geometry.z[1]),
            OFVector(geometry.x[0], geometry.y[6], geometry.z[0]),
            OFVector(geometry.x[0], geometry.y[6], geometry.z[1]),
            OFVector(geometry.x[0], geometry.yh[1], geometry.z[0]),
            OFVector(geometry.x[0], geometry.yh[1], geometry.z[1]),
            OFVector(geometry.x[0], geometry.y[11], geometry.z[0]),
            OFVector(geometry.x[0], geometry.y[11], geometry.z[1]),
            OFVector(geometry.x[0], geometry.y[12], geometry.z[0]),
            OFVector(geometry.x[0], geometry.y[12], geometry.z[1]),
            OFVector(geometry.x[1], geometry.y[0], geometry.z[0]),
            OFVector(geometry.x[1], geometry.y[0], geometry.z[1]),
            OFVector(geometry.x[1], geometry.y[1], geometry.z[0]),
            OFVector(geometry.x[1], geometry.y[1], geometry.z[1]),
            OFVector(geometry.x[1], geometry.yh[0], geometry.z[0]),
            OFVector(geometry.x[1], geometry.yh[0], geometry.z[1]),
            OFVector(geometry.x[1], geometry.y[6], geometry.z[0]),
            OFVector(geometry.x[1], geometry.y[6], geometry.z[1]),
            OFVector(geometry.x[1], geometry.yh[1], geometry.z[0]),
            OFVector(geometry.x[1], geometry.yh[1], geometry.z[1]),
            OFVector(geometry.x[1], geometry.y[11], geometry.z[0]),
            OFVector(geometry.x[1], geometry.y[11], geometry.z[1]),
            OFVector(geometry.x[1], geometry.y[12], geometry.z[0]),
            OFVector(geometry.x[1], geometry.y[12], geometry.z[1])]
        for i in range(geometry.n_wires):
            vertices.append(
                OFVector(
                    geometry.x[2+i*11], geometry.y[6], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[2+i*11], geometry.y[6], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[0], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[0], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[1], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[1], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[3], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[3], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[9], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[9], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[11], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[11], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[12], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[3+i*11], geometry.y[12], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[4+i*11], geometry.y[6], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[4+i*11], geometry.y[6], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[5+i*11], geometry.y[5], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[5+i*11], geometry.y[5], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[5+i*11], geometry.y[7], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[5+i*11], geometry.y[7], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[0], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[0], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[1], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[1], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[2], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[2], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[4], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[4], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[8], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[8], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[10], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[10], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[11], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[11], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[12], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[6+i*11], geometry.y[12], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[7+i*11], geometry.y[5], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[7+i*11], geometry.y[5], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[7+i*11], geometry.y[7], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[7+i*11], geometry.y[7], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[8+i*11], geometry.y[6], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[8+i*11], geometry.y[6], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[0], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[0], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[1], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[1], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[3], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[3], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[9], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[9], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[11], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[11], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[12], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[9+i*11], geometry.y[12], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[10+i*11], geometry.y[6], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[10+i*11], geometry.y[6], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[0], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[0], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[1], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[1], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.yh[0], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.yh[0], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[6], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[6], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.yh[1], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.yh[1], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[11], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[11], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[12], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[11+i*11], geometry.y[12], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[0], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[0], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[1], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[1], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.yh[0], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.yh[0], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[6], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[6], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.yh[1], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.yh[1], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[11], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[11], geometry.z[1]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[12], geometry.z[0]))
            vertices.append(
                OFVector(
                    geometry.x[12+i*11], geometry.y[12], geometry.z[1]))
        # outlet
        vertices_outlet = copy.deepcopy(vertices[-14:])
        for vertice in vertices_outlet:
            vertice[0] += geometry.x_outlet
        vertices += vertices_outlet
        # inlet
        vertices_inlet = copy.deepcopy(vertices[:14])
        for vertice in vertices_inlet:
            vertice[0] -= geometry.x_inlet
        vertices += vertices_inlet
        # wall
        vertices_wall = [
            OFVector(
                geometry.x[0] - geometry.x_inlet, 
                geometry.y[0] - geometry.d_wall/2.0, 
                geometry.z[0]),
            OFVector(
                geometry.x[0] - geometry.x_inlet, 
                geometry.y[0] - geometry.d_wall/2.0, 
                geometry.z[1]),
            OFVector(
                geometry.x[0] - geometry.x_inlet, 
                geometry.y[12] + geometry.d_wall/2.0, 
                geometry.z[0]),
            OFVector(
                geometry.x[0] - geometry.x_inlet, 
                geometry.y[12] + geometry.d_wall/2.0, 
                geometry.z[1]),
            OFVector(
                geometry.x[0], 
                geometry.y[0] - geometry.d_wall/2.0, 
                geometry.z[0]),
            OFVector(
                geometry.x[0], 
                geometry.y[0] - geometry.d_wall/2.0, 
                geometry.z[1]),
            OFVector(
                geometry.x[0], 
                geometry.y[12] + geometry.d_wall/2.0, 
                geometry.z[0]),
            OFVector(
                geometry.x[0], 
                geometry.y[12] + geometry.d_wall/2.0, 
                geometry.z[1]),
            OFVector(
                geometry.x[-1], 
                geometry.y[0] - geometry.d_wall/2.0, 
                geometry.z[0]),
            OFVector(
                geometry.x[-1], 
                geometry.y[0] - geometry.d_wall/2.0, 
                geometry.z[1]),
            OFVector(
                geometry.x[-1], 
                geometry.y[12] + geometry.d_wall/2.0, 
                geometry.z[0]),
            OFVector(
                geometry.x[-1], 
                geometry.y[12] + geometry.d_wall/2.0, 
                geometry.z[1]),
            OFVector(
                geometry.x[-1] + geometry.x_outlet, 
                geometry.y[0] - geometry.d_wall/2.0, 
                geometry.z[0]),
            OFVector(
                geometry.x[-1] + geometry.x_outlet, 
                geometry.y[0] - geometry.d_wall/2.0, 
                geometry.z[1]),
            OFVector(
                geometry.x[-1] + geometry.x_outlet, 
                geometry.y[12] + geometry.d_wall/2.0, 
                geometry.z[0]),
            OFVector(
                geometry.x[-1] + geometry.x_outlet, 
                geometry.y[12] + geometry.d_wall/2.0, 
                geometry.z[1])
        ]
        vertices += vertices_wall
        geometry.vertices = vertices
        # Create emptpy blockMeshDict dicts
        geometry.edges = {}
        geometry.blocks = {}
        geometry.boundary = {}
        geometry.mergePatchPairs = {}

    def calculateAirBlockMeshDict(self):
        """
        Calculates the blockMeshDict of the air regions
        """
        geometry = self.geometry
        # Create edges
        edges = []
        for i in range(geometry.n_wires):
            # 0'
            edges.append('arc')
            edges.append(42+i*84)
            edges.append(44+i*84)
            edges.append(
                [geometry.xp[1+i*8], geometry.yp[3], geometry.z[0]])
            edges.append('arc')
            edges.append(43+i*84)
            edges.append(45+i*84)
            edges.append(
                [geometry.xp[1+i*8], geometry.yp[3], geometry.z[1]])
            # 1'
            edges.append('arc')
            edges.append(44+i*84)
            edges.append(54+i*84)
            edges.append(
                [geometry.xp[3+i*8], geometry.yp[1], geometry.z[0]])
            edges.append('arc')
            edges.append(45+i*84)
            edges.append(55+i*84)
            edges.append(
                [geometry.xp[3+i*8], geometry.yp[1], geometry.z[1]])
            # 2'
            edges.append('arc')
            edges.append(54+i*84)
            edges.append(64+i*84)
            edges.append(
                [geometry.xp[4+i*8], geometry.yp[1], geometry.z[0]])
            edges.append('arc')
            edges.append(55+i*84)
            edges.append(65+i*84)
            edges.append(
                [geometry.xp[4+i*8], geometry.yp[1], geometry.z[1]])
            # 3'
            edges.append('arc')
            edges.append(64+i*84)
            edges.append(68+i*84)
            edges.append(
                [geometry.xp[6+i*8], geometry.yp[3], geometry.z[0]])
            edges.append('arc')
            edges.append(65+i*84)
            edges.append(69+i*84)
            edges.append(
                [geometry.xp[6+i*8], geometry.yp[3], geometry.z[1]])
            # 4'
            edges.append('arc')
            edges.append(68+i*84)
            edges.append(66+i*84)
            edges.append(
                [geometry.xp[6+i*8], geometry.yp[4], geometry.z[0]])
            edges.append('arc')
            edges.append(69+i*84)
            edges.append(67+i*84)
            edges.append(
                [geometry.xp[6+i*8], geometry.yp[4], geometry.z[1]])
            # 5'
            edges.append('arc')
            edges.append(66+i*84)
            edges.append(56+i*84)
            edges.append(
                [geometry.xp[4+i*8], geometry.yp[6], geometry.z[0]])
            edges.append('arc')
            edges.append(67+i*84)
            edges.append(57+i*84)
            edges.append(
                [geometry.xp[4+i*8], geometry.yp[6], geometry.z[1]])
            # 6'
            edges.append('arc')
            edges.append(56+i*84)
            edges.append(46+i*84)
            edges.append(
                [geometry.xp[3+i*8], geometry.yp[6], geometry.z[0]])
            edges.append('arc')
            edges.append(57+i*84)
            edges.append(47+i*84)
            edges.append(
                [geometry.xp[3+i*8], geometry.yp[6], geometry.z[1]])
            # 7'
            edges.append('arc')
            edges.append(46+i*84)
            edges.append(42+i*84)
            edges.append(
                [geometry.xp[1+i*8], geometry.yp[4], geometry.z[0]])
            edges.append('arc')
            edges.append(47+i*84)
            edges.append(43+i*84)
            edges.append(
                [geometry.xp[1+i*8], geometry.yp[4], geometry.z[1]])
            # 8'
            edges.append('arc')
            edges.append(28+i*84)
            edges.append(34+i*84)
            edges.append(
                [geometry.xp[0+i*8], geometry.yp[2], geometry.z[0]])
            edges.append('arc')
            edges.append(29+i*84)
            edges.append(35+i*84)
            edges.append(
                [geometry.xp[0+i*8], geometry.yp[2], geometry.z[1]])
            # 9'
            edges.append('arc')
            edges.append(34+i*84)
            edges.append(52+i*84)
            edges.append(
                [geometry.xp[2+i*8], geometry.yp[0], geometry.z[0]])
            edges.append('arc')
            edges.append(35+i*84)
            edges.append(53+i*84)
            edges.append(
                [geometry.xp[2+i*8], geometry.yp[0], geometry.z[1]])
            # 10'
            edges.append('arc')
            edges.append(52+i*84)
            edges.append(74+i*84)
            edges.append(
                [geometry.xp[5+i*8], geometry.yp[0], geometry.z[0]])
            edges.append('arc')
            edges.append(53+i*84)
            edges.append(75+i*84)
            edges.append(
                [geometry.xp[5+i*8], geometry.yp[0], geometry.z[1]])
            # 11'
            edges.append('arc')
            edges.append(74+i*84)
            edges.append(82+i*84)
            edges.append(
                [geometry.xp[7+i*8], geometry.yp[2], geometry.z[0]])
            edges.append('arc')
            edges.append(75+i*84)
            edges.append(83+i*84)
            edges.append(
                [geometry.xp[7+i*8], geometry.yp[2], geometry.z[1]])
            # 12'
            edges.append('arc')
            edges.append(82+i*84)
            edges.append(76+i*84)
            edges.append(
                [geometry.xp[7+i*8], geometry.yp[5], geometry.z[0]])
            edges.append('arc')
            edges.append(83+i*84)
            edges.append(77+i*84)
            edges.append(
                [geometry.xp[7+i*8], geometry.yp[5], geometry.z[1]])
            # 13'
            edges.append('arc')
            edges.append(76+i*84)
            edges.append(58+i*84)
            edges.append(
                [geometry.xp[5+i*8], geometry.yp[7], geometry.z[0]])
            edges.append('arc')
            edges.append(77+i*84)
            edges.append(59+i*84)
            edges.append(
                [geometry.xp[5+i*8], geometry.yp[7], geometry.z[1]])
            # 14'
            edges.append('arc')
            edges.append(58+i*84)
            edges.append(36+i*84)
            edges.append(
                [geometry.xp[2+i*8], geometry.yp[7], geometry.z[0]])
            edges.append('arc')
            edges.append(59+i*84)
            edges.append(37+i*84)
            edges.append(
                [geometry.xp[2+i*8], geometry.yp[7], geometry.z[1]])
            # 15'
            edges.append('arc')
            edges.append(36+i*84)
            edges.append(28+i*84)
            edges.append(
                [geometry.xp[0+i*8], geometry.yp[5], geometry.z[0]])
            edges.append('arc')
            edges.append(37+i*84)
            edges.append(29+i*84)
            edges.append(
                [geometry.xp[0+i*8], geometry.yp[5], geometry.z[1]])

        # Create blocks
        # 0 to 3
        blocks = [
            'hex',
            [2, 16, 18, 4, 3, 17, 19, 5],
            [geometry.nx[0], geometry.ny1, geometry.nz],
            'simpleGrading',
            '(1 1 1)',
            'hex',
            [4, 18, 20, 6, 5, 19, 21, 7],
            [geometry.nx[0], geometry.ny1, geometry.nz],
            'simpleGrading',
            '(1 1 1)',
            'hex',
            [6, 20, 22, 8, 7, 21, 23, 9],
            [geometry.nx[0], geometry.ny1, geometry.nz],
            'simpleGrading',
            '(1 1 1)',
            'hex',
            [8, 22, 24, 10, 9, 23, 25, 11],
            [geometry.nx[0], geometry.ny1, geometry.nz],
            'simpleGrading',
            '(1 1 1)']
        for i in range(geometry.n_wires):
            # 4
            blocks.append('hex')
            blocks.append([
                16+i*84, 32+i*84, 34+i*84, 18+i*84, 
                17+i*84, 33+i*84, 35+i*84, 19+i*84])
            blocks.append([geometry.nx[1], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 5
            blocks.append('hex')
            blocks.append([
                18+i*84, 34+i*84, 28+i*84, 20+i*84,
                19+i*84, 35+i*84, 29+i*84, 21+i*84])
            blocks.append([geometry.nx[1], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 6
            blocks.append('hex')
            blocks.append([
                20+i*84, 28+i*84, 36+i*84, 22+i*84,
                21+i*84, 29+i*84, 37+i*84, 23+i*84])
            blocks.append([geometry.nx[1], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 7
            blocks.append('hex')
            blocks.append([
                22+i*84, 36+i*84, 38+i*84, 24+i*84,
                23+i*84, 37+i*84, 39+i*84, 25+i*84])
            blocks.append([geometry.nx[1], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 8
            blocks.append('hex')
            blocks.append([
                32+i*84, 50+i*84, 52+i*84, 34+i*84,
                33+i*84, 51+i*84, 53+i*84, 35+i*84])
            blocks.append([geometry.nx[2], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 9
            blocks.append('hex')
            blocks.append([
                34+i*84, 52+i*84, 54+i*84, 44+i*84,
                35+i*84, 53+i*84, 55+i*84, 45+i*84])
            blocks.append([geometry.nx[2], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 10
            blocks.append('hex')
            blocks.append([
                44+i*84, 42+i*84, 28+i*84, 34+i*84,
                45+i*84, 43+i*84, 29+i*84, 35+i*84])
            blocks.append([geometry.nx[2], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 11
            blocks.append('hex')
            blocks.append([
                42+i*84, 46+i*84, 36+i*84, 28+i*84,
                43+i*84, 47+i*84, 37+i*84, 29+i*84])
            blocks.append([geometry.nx[2], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 12
            blocks.append('hex')
            blocks.append([
                46+i*84, 56+i*84, 58+i*84, 36+i*84,
                47+i*84, 57+i*84, 59+i*84, 37+i*84])
            blocks.append([geometry.nx[2], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 13
            blocks.append('hex')
            blocks.append([
                36+i*84, 58+i*84, 60+i*84, 38+i*84,
                37+i*84, 59+i*84, 61+i*84, 39+i*84])
            blocks.append([geometry.nx[2], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 14
            blocks.append('hex')
            blocks.append([
                50+i*84, 72+i*84, 74+i*84, 52+i*84,
                51+i*84, 73+i*84, 75+i*84, 53+i*84])
            blocks.append([geometry.nx[3], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 15
            blocks.append('hex')
            blocks.append([
                52+i*84, 74+i*84, 64+i*84, 54+i*84,
                53+i*84, 75+i*84, 65+i*84, 55+i*84])
            blocks.append([geometry.nx[3], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 16
            blocks.append('hex')
            blocks.append([
                74+i*84, 82+i*84, 68+i*84, 64+i*84,
                75+i*84, 83+i*84, 69+i*84, 65+i*84])
            blocks.append([geometry.nx[3], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 17
            blocks.append('hex')
            blocks.append([
                82+i*84, 76+i*84, 66+i*84, 68+i*84,
                83+i*84, 77+i*84, 67+i*84, 69+i*84])
            blocks.append([geometry.nx[3], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 18
            blocks.append('hex')
            blocks.append([
                76+i*84, 58+i*84, 56+i*84, 66+i*84,
                77+i*84, 59+i*84, 57+i*84, 67+i*84])
            blocks.append([geometry.nx[3], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 19
            blocks.append('hex')
            blocks.append([
                58+i*84, 76+i*84, 78+i*84, 60+i*84,
                59+i*84, 77+i*84, 79+i*84, 61+i*84])
            blocks.append([geometry.nx[3], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 20
            blocks.append('hex')
            blocks.append([
                72+i*84, 86+i*84, 88+i*84, 74+i*84,
                73+i*84, 87+i*84, 89+i*84, 75+i*84])
            blocks.append([geometry.nx[4], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 21
            blocks.append('hex')
            blocks.append([
                74+i*84, 88+i*84, 90+i*84, 82+i*84,
                75+i*84, 89+i*84, 91+i*84, 83+i*84])
            blocks.append([geometry.nx[4], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 22
            blocks.append('hex')
            blocks.append([
                82+i*84, 90+i*84, 92+i*84, 76+i*84,
                83+i*84, 91+i*84, 93+i*84, 77+i*84])
            blocks.append([geometry.nx[4], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 23
            blocks.append('hex')
            blocks.append([
                76+i*84, 92+i*84, 94+i*84, 78+i*84,
                77+i*84, 93+i*84, 95+i*84, 79+i*84])
            blocks.append([geometry.nx[4], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 24
            blocks.append('hex')
            blocks.append([
                86+i*84, 100+i*84, 102+i*84, 88+i*84,
                87+i*84, 101+i*84, 103+i*84, 89+i*84])
            blocks.append([geometry.nx[5+i], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 25
            blocks.append('hex')
            blocks.append([
                88+i*84, 102+i*84, 104+i*84, 90+i*84,
                89+i*84, 103+i*84, 105+i*84, 91+i*84])
            blocks.append([geometry.nx[5+i], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 26
            blocks.append('hex')
            blocks.append([
                90+i*84, 104+i*84, 106+i*84, 92+i*84,
                91+i*84, 105+i*84, 107+i*84, 93+i*84])
            blocks.append([geometry.nx[5+i], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
            # 27
            blocks.append('hex')
            blocks.append([
                92+i*84, 106+i*84, 108+i*84, 94+i*84,
                93+i*84, 107+i*84, 109+i*84, 95+i*84])
            blocks.append([geometry.nx[5+i], geometry.ny1, geometry.nz])
            blocks.append('simpleGrading')
            blocks.append('(1 1 1)')
        # outlet
        num_vertices = 28 + geometry.n_wires * 84
        blocks_outlet = copy.deepcopy(blocks[-20:])
        for i in range(1,20,5):
            for ii in range(8):
                blocks_outlet[i][ii] += 14
        for i in range(2, 20, 5):
            blocks_outlet[i][0] = geometry.nx_outlet
        if geometry.grading_inlet_outlet:
            for i in range(4, 20, 5):
                blocks_outlet[i] = (f'({geometry.grading_outlet} 1 1)')
        # solution blocks
        blocks_outlet.append('hex')
        blocks_outlet.append([
            num_vertices - 14,
            num_vertices + 0,
            num_vertices + 2,
            num_vertices - 12,
            num_vertices - 13,
            num_vertices + 1,
            num_vertices + 3,
            num_vertices - 11])
        blocks_outlet.append(
            [geometry.nx_outlet, geometry.ny_sol_in_out, geometry.nz])
        blocks_outlet.append('simpleGrading')
        if geometry.grading_inlet_outlet:
            blocks_outlet.append(f'({geometry.grading_outlet} 1 1)')
        else:
            blocks_outlet.append('(1 1 1)')
        blocks_outlet.append('hex')
        blocks_outlet.append([
            num_vertices - 4,
            num_vertices + 10,
            num_vertices + 12,
            num_vertices - 2,
            num_vertices - 3,
            num_vertices + 11,
            num_vertices + 13,
            num_vertices - 1])
        blocks_outlet.append(
            [geometry.nx_outlet, geometry.ny_sol_in_out, geometry.nz])
        blocks_outlet.append('simpleGrading')
        if geometry.grading_inlet_outlet:
            blocks_outlet.append(f'({geometry.grading_outlet} 1 1)')
        else:
            blocks_outlet.append('(1 1 1)')
        # wall blocks
        blocks_outlet.append('hex')
        blocks_outlet.append([
            num_vertices + 36,
            num_vertices + 40,
            num_vertices + 0,
            num_vertices - 14,
            num_vertices + 37,
            num_vertices + 41,
            num_vertices + 1,
            num_vertices - 13])
        blocks_outlet.append(
            [geometry.nx_outlet, geometry.ny_wall, geometry.nz])
        blocks_outlet.append('simpleGrading')
        if geometry.grading_inlet_outlet:
            blocks_outlet.append(f'({geometry.grading_outlet} 1 1)')
        else:
            blocks_outlet.append('(1 1 1)')
        blocks_outlet.append('hex')
        blocks_outlet.append([
            num_vertices - 2,
            num_vertices + 12,
            num_vertices + 42,
            num_vertices + 38,
            num_vertices - 1,
            num_vertices + 13,
            num_vertices + 43,
            num_vertices + 39])
        blocks_outlet.append(
            [geometry.nx_outlet, geometry.ny_wall, geometry.nz])
        blocks_outlet.append('simpleGrading')
        if geometry.grading_inlet_outlet:
            blocks_outlet.append(f'({geometry.grading_outlet} 1 1)')
        else:
            blocks_outlet.append('(1 1 1)')
        blocks += blocks_outlet
        # inlet
        blocks_inlet = copy.deepcopy(blocks[:20])
        for i in range(1,20,5):
            for ii in range(8):
                if ii in [1, 2, 5, 6]:
                    blocks_inlet[i][ii] -= 14
                else:
                    blocks_inlet[i][ii] += (num_vertices + 14)
        for i in range(2, 20, 5):
            blocks_inlet[i][0] = geometry.nx_inlet
        if geometry.grading_inlet_outlet:
            for i in range(4, 20, 5):
                blocks_inlet[i] = (f'({geometry.grading_inlet} 1 1)')
        # solution blocks
        blocks_inlet.append('hex')
        blocks_inlet.append([
            num_vertices + 14,
            0,
            2,
            num_vertices + 16,
            num_vertices + 15,
            1,
            3,
            num_vertices + 17])
        blocks_inlet.append(
            [geometry.nx_inlet, geometry.ny_sol_in_out, geometry.nz])
        blocks_inlet.append('simpleGrading')
        if geometry.grading_inlet_outlet:
            blocks_inlet.append(f'({geometry.grading_inlet} 1 1)')
        else:
            blocks_inlet.append('(1 1 1)')
        blocks_inlet.append('hex')
        blocks_inlet.append([
            num_vertices + 24,
            10,
            12,
            num_vertices + 26,
            num_vertices + 25,
            11,
            13,
            num_vertices + 27])
        blocks_inlet.append(
            [geometry.nx_inlet, geometry.ny_sol_in_out, geometry.nz])
        blocks_inlet.append('simpleGrading')
        if geometry.grading_inlet_outlet:
            blocks_inlet.append(f'({geometry.grading_inlet} 1 1)')
        else:
            blocks_inlet.append('(1 1 1)')
        # wall blocks
        blocks_inlet.append('hex')
        blocks_inlet.append([
            num_vertices + 28,
            num_vertices + 32,
            0,
            num_vertices + 14,
            num_vertices + 29,
            num_vertices + 33,
            1,
            num_vertices + 15])
        blocks_inlet.append(
            [geometry.nx_inlet, geometry.ny_wall, geometry.nz])
        blocks_inlet.append('simpleGrading')
        if geometry.grading_inlet_outlet:
            blocks_inlet.append(f'({geometry.grading_inlet} 1 1)')
        else:
            blocks_inlet.append('(1 1 1)')
        blocks_inlet.append('hex')
        blocks_inlet.append([
            num_vertices + 26,
            12,
            num_vertices + 34,
            num_vertices + 30,
            num_vertices + 27,
            13,
            num_vertices + 35,
            num_vertices + 31])
        blocks_inlet.append(
            [geometry.nx_inlet, geometry.ny_wall, geometry.nz])
        blocks_inlet.append('simpleGrading')
        if geometry.grading_inlet_outlet:
            blocks_inlet.append(f'({geometry.grading_inlet} 1 1)')
        else:
            blocks_inlet.append('(1 1 1)')
        blocks += blocks_inlet
        # Create boundary
        # Define faces
        inlet_faces = [
                [2, 4, 5, 3],
                [4, 6, 7, 5],
                [6, 8, 9, 7],
                [8, 10, 11, 9]]
        inlet_faces = [
            list(map(lambda x: x + num_vertices + 14, face)) 
            for face in inlet_faces]
        # solution
        inlet_faces.append([
            num_vertices + 14,
            num_vertices + 16,
            num_vertices + 17,
            num_vertices + 15])
        inlet_faces.append([
            num_vertices + 24,
            num_vertices + 26,
            num_vertices + 27,
            num_vertices + 25])
        # wall
        inlet_faces.append([
            num_vertices + 28,
            num_vertices + 14,
            num_vertices + 15,
            num_vertices + 29])
        inlet_faces.append([
            num_vertices + 26,
            num_vertices + 30,
            num_vertices + 31,
            num_vertices + 27])
        outlet_faces = [
                [16, 18, 19, 17],
                [18, 20, 21, 19],
                [20, 22, 23, 21],
                [22, 24, 25, 23]]
        # Add geometry.n_wires*84 to all vortices of outlet_faces
        outlet_faces = [
            list(map(lambda x: x + geometry.n_wires*84 + 14, face)) 
            for face in outlet_faces]
        # solution
        outlet_faces.append([
            num_vertices + 0,
            num_vertices + 2,
            num_vertices + 3,
            num_vertices + 1])
        outlet_faces.append([
            num_vertices + 10,
            num_vertices + 12,
            num_vertices + 13,
            num_vertices + 11])
        # wall
        outlet_faces.append([
            num_vertices + 40,
            num_vertices + 0,
            num_vertices + 1,
            num_vertices + 41])
        outlet_faces.append([
            num_vertices + 12,
            num_vertices + 42,
            num_vertices + 43,
            num_vertices + 13])
        air_to_solution1_faces = [
            [2, 16, 17, 3]]
        for i in range(geometry.n_wires):
            air_to_solution1_faces.append([
                16+i*84, 32+i*84, 33+i*84, 17+i*84])
            air_to_solution1_faces.append([
                32+i*84, 50+i*84, 51+i*84, 33+i*84])
            air_to_solution1_faces.append([
                50+i*84, 72+i*84, 73+i*84, 51+i*84])
            air_to_solution1_faces.append([
                72+i*84, 86+i*84, 87+i*84, 73+i*84])
            air_to_solution1_faces.append([
                86+i*84, 100+i*84, 101+i*84, 87+i*84])
        air_to_solution2_faces = [
            [10, 24, 25, 11]]
        for i in range(geometry.n_wires):
            air_to_solution2_faces.append([
                24+i*84, 38+i*84, 39+i*84, 25+i*84])
            air_to_solution2_faces.append([
                38+i*84, 60+i*84, 61+i*84, 39+i*84])
            air_to_solution2_faces.append([
                60+i*84, 78+i*84, 79+i*84, 61+i*84])
            air_to_solution2_faces.append([
                78+i*84, 94+i*84, 95+i*84, 79+i*84])
            air_to_solution2_faces.append([
                94+i*84, 108+i*84, 109+i*84, 95+i*84])
        wire_faces = []
        for i in range(geometry.n_wires):
            wire_faces.append(
                [42+i*84, 44+i*84, 45+i*84, 43+i*84])
            wire_faces.append(
                [44+i*84, 54+i*84, 55+i*84, 45+i*84])
            wire_faces.append(
                [54+i*84, 64+i*84, 65+i*84, 55+i*84])
            wire_faces.append(
                [64+i*84, 68+i*84, 69+i*84, 65+i*84])
            wire_faces.append(
                [68+i*84, 66+i*84, 67+i*84, 69+i*84])
            wire_faces.append(
                [66+i*84, 56+i*84, 57+i*84, 67+i*84])
            wire_faces.append(
                [56+i*84, 46+i*84, 47+i*84, 57+i*84])
            wire_faces.append(
                [46+i*84, 42+i*84, 43+i*84, 47+i*84])
        top_bottom_faces = [
            [2, 16, 18, 4],
            [3, 17, 19, 5],
            [4, 18, 20, 6],
            [5, 19, 21, 7],
            [6, 20, 22, 8],
            [7, 21, 23, 9],
            [8, 22, 24, 10],
            [9, 23, 25, 11]]
        for i in range(geometry.n_wires):
            top_bottom_faces.append(
                [16+i*84, 32+i*84, 34+i*84, 18+i*84])
            top_bottom_faces.append(
                [17+i*84, 33+i*84, 35+i*84, 19+i*84])
            top_bottom_faces.append(
                [18+i*84, 34+i*84, 28+i*84, 20+i*84])
            top_bottom_faces.append(
                [19+i*84, 35+i*84, 29+i*84, 21+i*84])
            top_bottom_faces.append(
                [20+i*84, 28+i*84, 36+i*84, 22+i*84])
            top_bottom_faces.append(
                [21+i*84, 29+i*84, 37+i*84, 23+i*84])
            top_bottom_faces.append(
                [22+i*84, 36+i*84, 38+i*84, 24+i*84])
            top_bottom_faces.append(
                [23+i*84, 37+i*84, 39+i*84, 25+i*84])
            top_bottom_faces.append(
                [32+i*84, 50+i*84, 52+i*84, 34+i*84])
            top_bottom_faces.append(
                [33+i*84, 51+i*84, 53+i*84, 35+i*84])
            top_bottom_faces.append(
                [50+i*84, 72+i*84, 74+i*84, 52+i*84])
            top_bottom_faces.append(
                [51+i*84, 73+i*84, 75+i*84, 53+i*84])
            top_bottom_faces.append(
                [28+i*84, 34+i*84, 44+i*84, 42+i*84])
            top_bottom_faces.append(
                [29+i*84, 35+i*84, 45+i*84, 43+i*84])
            top_bottom_faces.append(
                [34+i*84, 52+i*84, 54+i*84, 44+i*84])
            top_bottom_faces.append(
                [35+i*84, 53+i*84, 55+i*84, 45+i*84])
            top_bottom_faces.append(
                [52+i*84, 74+i*84, 64+i*84, 54+i*84])
            top_bottom_faces.append(
                [53+i*84, 75+i*84, 65+i*84, 55+i*84])
            top_bottom_faces.append(
                [74+i*84, 82+i*84, 68+i*84, 64+i*84])
            top_bottom_faces.append(
                [75+i*84, 83+i*84, 69+i*84, 65+i*84])
            top_bottom_faces.append(
                [82+i*84, 76+i*84, 66+i*84, 68+i*84])
            top_bottom_faces.append(
                [83+i*84, 77+i*84, 67+i*84, 69+i*84])
            top_bottom_faces.append(
                [76+i*84, 58+i*84, 56+i*84, 66+i*84])
            top_bottom_faces.append(
                [77+i*84, 59+i*84, 57+i*84, 67+i*84])
            top_bottom_faces.append(
                [58+i*84, 36+i*84, 46+i*84, 56+i*84])
            top_bottom_faces.append(
                [59+i*84, 37+i*84, 47+i*84, 57+i*84])
            top_bottom_faces.append(
                [36+i*84, 28+i*84, 42+i*84, 46+i*84])
            top_bottom_faces.append(
                [37+i*84, 29+i*84, 43+i*84, 47+i*84])
            top_bottom_faces.append(
                [36+i*84, 58+i*84, 60+i*84, 38+i*84])
            top_bottom_faces.append(
                [37+i*84, 59+i*84, 61+i*84, 39+i*84])
            top_bottom_faces.append(
                [58+i*84, 76+i*84, 78+i*84, 60+i*84])
            top_bottom_faces.append(
                [59+i*84, 77+i*84, 79+i*84, 61+i*84])
            top_bottom_faces.append(
                [72+i*84, 86+i*84, 88+i*84, 74+i*84])
            top_bottom_faces.append(
                [73+i*84, 87+i*84, 89+i*84, 75+i*84])
            top_bottom_faces.append(
                [74+i*84, 88+i*84, 90+i*84, 82+i*84])
            top_bottom_faces.append(
                [75+i*84, 89+i*84, 91+i*84, 83+i*84])
            top_bottom_faces.append(
                [82+i*84, 90+i*84, 92+i*84, 76+i*84])
            top_bottom_faces.append(
                [83+i*84, 91+i*84, 93+i*84, 77+i*84])
            top_bottom_faces.append(
                [76+i*84, 92+i*84, 94+i*84, 78+i*84])
            top_bottom_faces.append(
                [77+i*84, 93+i*84, 95+i*84, 79+i*84])
            top_bottom_faces.append(
                [86+i*84, 100+i*84, 102+i*84, 88+i*84])
            top_bottom_faces.append(
                [87+i*84, 101+i*84, 103+i*84, 89+i*84])
            top_bottom_faces.append(
                [88+i*84, 102+i*84, 104+i*84, 90+i*84])
            top_bottom_faces.append(
                [89+i*84, 103+i*84, 105+i*84, 91+i*84])
            top_bottom_faces.append(
                [90+i*84, 104+i*84, 106+i*84, 92+i*84])
            top_bottom_faces.append(
                [91+i*84, 105+i*84, 107+i*84, 93+i*84])
            top_bottom_faces.append(
                [92+i*84, 106+i*84, 108+i*84, 94+i*84])
            top_bottom_faces.append(
                [93+i*84, 107+i*84, 109+i*84, 95+i*84])
        # outlet
        top_bottom_faces.append([
            num_vertices + 36,
            num_vertices + 40,
            num_vertices + 0,
            num_vertices - 14])
        top_bottom_faces.append([
            num_vertices + 37,
            num_vertices + 41,
            num_vertices + 1,
            num_vertices - 13])
        top_bottom_faces.append([
            num_vertices - 14,
            num_vertices + 0,
            num_vertices + 2,
            num_vertices - 12])
        top_bottom_faces.append([
            num_vertices - 13,
            num_vertices + 1,
            num_vertices + 3,
            num_vertices - 11])
        top_bottom_faces.append([
            num_vertices - 12,
            num_vertices + 2,
            num_vertices + 4,
            num_vertices - 10])
        top_bottom_faces.append([
            num_vertices - 11,
            num_vertices + 3,
            num_vertices + 5,
            num_vertices - 9])
        top_bottom_faces.append([
            num_vertices - 10,
            num_vertices + 4,
            num_vertices + 6,
            num_vertices - 8])
        top_bottom_faces.append([
            num_vertices - 9,
            num_vertices + 5,
            num_vertices + 7,
            num_vertices - 7])
        top_bottom_faces.append([
            num_vertices - 8,
            num_vertices + 6,
            num_vertices + 8,
            num_vertices - 6])
        top_bottom_faces.append([
            num_vertices - 7,
            num_vertices + 7,
            num_vertices + 9,
            num_vertices - 5])
        top_bottom_faces.append([
            num_vertices - 6,
            num_vertices + 8,
            num_vertices + 10,
            num_vertices - 4])
        top_bottom_faces.append([
            num_vertices - 5,
            num_vertices + 9,
            num_vertices + 11,
            num_vertices - 3])
        top_bottom_faces.append([
            num_vertices - 4,
            num_vertices + 10,
            num_vertices + 12,
            num_vertices - 2])
        top_bottom_faces.append([
            num_vertices - 3,
            num_vertices + 11,
            num_vertices + 13,
            num_vertices - 1])
        top_bottom_faces.append([
            num_vertices - 2,
            num_vertices + 12,
            num_vertices + 42,
            num_vertices + 38])
        top_bottom_faces.append([
            num_vertices - 1,
            num_vertices + 13,
            num_vertices + 43,
            num_vertices + 39])
        # inlet
        top_bottom_faces.append([
            num_vertices + 28,
            num_vertices + 32,
            0,
            num_vertices + 14])
        top_bottom_faces.append([
            num_vertices + 29,
            num_vertices + 33,
            1,
            num_vertices + 15])
        top_bottom_faces.append([
            num_vertices + 14,
            0,
            2,
            num_vertices + 16])
        top_bottom_faces.append([
            num_vertices + 15,
            1,
            3,
            num_vertices + 17])
        top_bottom_faces.append([
            num_vertices + 16,
            2,
            4,
            num_vertices + 18])
        top_bottom_faces.append([
            num_vertices + 17,
            3,
            5,
            num_vertices + 19])
        top_bottom_faces.append([
            num_vertices + 18,
            4,
            6,
            num_vertices + 20])
        top_bottom_faces.append([
            num_vertices + 19,
            5,
            7,
            num_vertices + 21])
        top_bottom_faces.append([
            num_vertices + 20,
            6,
            8,
            num_vertices + 22])
        top_bottom_faces.append([
            num_vertices + 21,
            7,
            9,
            num_vertices + 23])
        top_bottom_faces.append([
            num_vertices + 22,
            8,
            10,
            num_vertices + 24])
        top_bottom_faces.append([
            num_vertices + 23,
            9,
            11,
            num_vertices + 25])
        top_bottom_faces.append([
            num_vertices + 24,
            10,
            12,
            num_vertices + 26])
        top_bottom_faces.append([
            num_vertices + 25,
            11,
            13,
            num_vertices + 27])
        top_bottom_faces.append([
            num_vertices + 26,
            12,
            num_vertices + 34,
            num_vertices + 30])
        top_bottom_faces.append([
            num_vertices + 27,
            13,
            num_vertices + 35,
            num_vertices + 31])
        # wall
        wall_faces = []
        # inlet
        wall_faces.append([
            num_vertices + 32,
            0,
            1,
            num_vertices + 33])
        wall_faces.append([
            0,
            2,
            3,
            1])
        wall_faces.append([
            10,
            12,
            13,
            11])
        wall_faces.append([
            12,
            num_vertices + 34,
            num_vertices + 35,
            13])
        # outlet
        wall_faces.append([
            num_vertices + 36,
            num_vertices - 14,
            num_vertices - 13,
            num_vertices + 37])
        wall_faces.append([
            num_vertices - 14,
            num_vertices - 12,
            num_vertices - 11,
            num_vertices - 13])
        wall_faces.append([
            num_vertices - 4,
            num_vertices - 2,
            num_vertices - 1,
            num_vertices - 3])
        wall_faces.append([
            num_vertices - 2,
            num_vertices + 38,
            num_vertices + 39,
            num_vertices - 1])
        # symmetry
        symmetry_faces = []
        # inlet
        symmetry_faces.append([
            num_vertices + 28,
            num_vertices + 32,
            num_vertices + 33,
            num_vertices + 29])
        symmetry_faces.append([
            num_vertices + 30,
            num_vertices + 34,
            num_vertices + 35,
            num_vertices + 31])
        # outlet
        symmetry_faces.append([
            num_vertices + 36,
            num_vertices + 40,
            num_vertices + 41,
            num_vertices + 37])
        symmetry_faces.append([
            num_vertices + 38,
            num_vertices + 42,
            num_vertices + 43,
            num_vertices + 39])
        boundary = [
            'inlet',
            {'type': 'patch', 'faces': inlet_faces},
            'outlet',
            {'type': 'patch', 'faces': outlet_faces},
            'air_to_solution1',
            {'type': 'patch',
            'faces': air_to_solution1_faces},
            'air_to_solution2',
            {'type': 'patch',
            'faces': air_to_solution2_faces},
            'wire',
            {'type': 'patch',
            'faces': wire_faces},
            'top_bottom',
            {'type': 'patch',
            'faces': top_bottom_faces},
            'wall',
            {'type': 'wall',
             'faces': wall_faces},
            'symmetry',
            {'type': 'symmetry',
             'faces': symmetry_faces}]
       
        # Create mergePatchPairs
        mergePatchPairs = []

        for region in self.solver.airRegions:
            geometry.edges[region] = edges
            geometry.blocks[region] = blocks
            geometry.boundary[region] = boundary
            geometry.mergePatchPairs[region] = mergePatchPairs

    def calculateSolutionBlockMeshDict(self):
        """
        Calculates the blockMeshDict of the air region
        """
        geometry = self.geometry
        # Create edges
        edges = []
       
        # Create mergePatchPairs
        mergePatchPairs = []

        for region in self.solver.solutionRegions:
            if region == 'solution1':
                # Create blocks
                # 0
                blocks = [
                    'hex',
                    [0, 14, 16, 2, 1, 15, 17, 3],
                    [geometry.nx[0], geometry.ny_sol, geometry.nz],
                    'simpleGrading',
                    f'(1 {1/self.geometry.solution.grading} 1)']
                for i in range(geometry.n_wires):
                    # 1
                    blocks.append('hex')
                    blocks.append([
                        14+i*84, 30+i*84, 32+i*84, 16+i*84, 
                        15+i*84, 31+i*84, 33+i*84, 17+i*84])
                    blocks.append(
                        [geometry.nx[1], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {1/self.geometry.solution.grading} 1)')
                    # 2
                    blocks.append('hex')
                    blocks.append([
                        30+i*84, 48+i*84, 50+i*84, 32+i*84,
                        31+i*84, 49+i*84, 51+i*84, 33+i*84])
                    blocks.append(
                        [geometry.nx[2], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {1/self.geometry.solution.grading} 1)')
                    # 3
                    blocks.append('hex')
                    blocks.append([
                        48+i*84, 70+i*84, 72+i*84, 50+i*84,
                        49+i*84, 71+i*84, 73+i*84, 51+i*84])
                    blocks.append(
                        [geometry.nx[3], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {1/self.geometry.solution.grading} 1)')
                    # 4
                    blocks.append('hex')
                    blocks.append([
                        70+i*84, 84+i*84, 86+i*84, 72+i*84,
                        71+i*84, 85+i*84, 87+i*84, 73+i*84])
                    blocks.append(
                        [geometry.nx[4], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {1/self.geometry.solution.grading} 1)')
                    # 5
                    blocks.append('hex')
                    blocks.append([
                        84+i*84, 98+i*84, 100+i*84, 86+i*84,
                        85+i*84, 99+i*84, 101+i*84, 87+i*84])
                    blocks.append(
                        [geometry.nx[5+i], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {1/self.geometry.solution.grading} 1)')
                
                # Create boundary
                # Define faces
                inlet_faces = [[0, 14, 16, 2]]
                for i in range(geometry.n_wires):
                    inlet_faces.append([14+i*84, 30+i*84, 32+i*84, 16+i*84])
                    inlet_faces.append([30+i*84, 48+i*84, 50+i*84, 32+i*84])
                    inlet_faces.append([48+i*84, 70+i*84, 72+i*84, 50+i*84])
                    inlet_faces.append([70+i*84, 84+i*84, 86+i*84, 72+i*84])
                    inlet_faces.append([84+i*84, 98+i*84, 100+i*84, 86+i*84])
                outlet_faces = [[1, 15, 17, 3]]
                for i in range(geometry.n_wires):
                    outlet_faces.append([15+i*84, 31+i*84, 33+i*84, 17+i*84])
                    outlet_faces.append([31+i*84, 49+i*84, 51+i*84, 33+i*84])
                    outlet_faces.append([49+i*84, 71+i*84, 73+i*84, 51+i*84])
                    outlet_faces.append([71+i*84, 85+i*84, 87+i*84, 73+i*84])
                    outlet_faces.append([85+i*84, 99+i*84, 101+i*84, 87+i*84])
                solution_to_air_faces = [[2, 16, 17, 3]]
                for i in range(geometry.n_wires):
                    solution_to_air_faces.append([
                        16+i*84, 32+i*84, 33+i*84, 17+i*84])
                    solution_to_air_faces.append([
                        32+i*84, 50+i*84, 51+i*84, 33+i*84])
                    solution_to_air_faces.append([
                        50+i*84, 72+i*84, 73+i*84, 51+i*84])
                    solution_to_air_faces.append([
                        72+i*84, 86+i*84, 87+i*84, 73+i*84])
                    solution_to_air_faces.append([
                        86+i*84, 100+i*84, 101+i*84, 87+i*84])
                wall_faces = [[0, 14, 15, 1]]
                for i in range(geometry.n_wires):
                    wall_faces.append([14+i*84, 30+i*84, 31+i*84, 15+i*84])
                    wall_faces.append([30+i*84, 48+i*84, 49+i*84, 31+i*84])
                    wall_faces.append([48+i*84, 70+i*84, 71+i*84, 49+i*84])
                    wall_faces.append([70+i*84, 84+i*84, 85+i*84, 71+i*84])
                    wall_faces.append([84+i*84, 98+i*84, 99+i*84, 85+i*84])
                front_back_faces = [
                    [0, 2, 3, 1], 
                        [14+geometry.n_wires*84, 
                        16+geometry.n_wires*84, 
                        17+geometry.n_wires*84, 
                        15+geometry.n_wires*84]]

                boundary = [
                    'inlet',
                    {'type': 'patch', 'faces': inlet_faces},
                    'outlet',
                    {'type': 'patch', 'faces': outlet_faces},
                    'wall',
                    {'type': 'patch',
                    'faces': wall_faces},
                    'front_back',
                    {'type': 'patch',
                    'faces': front_back_faces},
                    'solution_to_air',
                    {'type': 'patch',
                    'faces': solution_to_air_faces}]
            elif region == 'solution2':
                # Create blocks
                # 0
                blocks = [
                    'hex',
                    [10, 24, 26, 12, 11, 25, 27, 13],
                    [geometry.nx[0], geometry.ny_sol, geometry.nz],
                    'simpleGrading',
                    f'(1 {self.geometry.solution.grading} 1)']
                for i in range(geometry.n_wires):
                    # 1
                    blocks.append('hex')
                    blocks.append([
                        24+i*84, 38+i*84, 40+i*84, 26+i*84, 
                        25+i*84, 39+i*84, 41+i*84, 27+i*84])
                    blocks.append(
                        [geometry.nx[1], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {self.geometry.solution.grading} 1)')
                    # 2
                    blocks.append('hex')
                    blocks.append([
                        38+i*84, 60+i*84, 62+i*84, 40+i*84,
                        39+i*84, 61+i*84, 63+i*84, 41+i*84])
                    blocks.append(
                        [geometry.nx[2], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {self.geometry.solution.grading} 1)')
                    # 3
                    blocks.append('hex')
                    blocks.append([
                        60+i*84, 78+i*84, 80+i*84, 62+i*84,
                        61+i*84, 79+i*84, 81+i*84, 63+i*84])
                    blocks.append(
                        [geometry.nx[3], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {self.geometry.solution.grading} 1)')
                    # 4
                    blocks.append('hex')
                    blocks.append([
                        78+i*84, 94+i*84, 96+i*84, 80+i*84,
                        79+i*84, 95+i*84, 97+i*84, 81+i*84])
                    blocks.append(
                        [geometry.nx[4], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {self.geometry.solution.grading} 1)')
                    # 5
                    blocks.append('hex')
                    blocks.append([
                        94+i*84, 108+i*84, 110+i*84, 96+i*84,
                        95+i*84, 109+i*84, 111+i*84, 97+i*84])
                    blocks.append(
                        [geometry.nx[5+i], geometry.ny_sol, geometry.nz])
                    blocks.append('simpleGrading')
                    blocks.append(f'(1 {self.geometry.solution.grading} 1)')

                # Create boundary
                # Define faces
                inlet_faces = [[10, 12, 24, 26]]
                for i in range(geometry.n_wires):
                    inlet_faces.append([24+i*84, 38+i*84, 40+i*84, 26+i*84])
                    inlet_faces.append([38+i*84, 60+i*84, 62+i*84, 40+i*84])
                    inlet_faces.append([60+i*84, 78+i*84, 80+i*84, 62+i*84])
                    inlet_faces.append([78+i*84, 94+i*84, 96+i*84, 80+i*84])
                    inlet_faces.append([94+i*84, 108+i*84, 110+i*84, 96+i*84])
                outlet_faces = [[11, 13, 25, 27]]
                for i in range(geometry.n_wires):
                    outlet_faces.append([25+i*84, 39+i*84, 41+i*84, 27+i*84])
                    outlet_faces.append([39+i*84, 61+i*84, 63+i*84, 41+i*84])
                    outlet_faces.append([61+i*84, 79+i*84, 81+i*84, 63+i*84])
                    outlet_faces.append([79+i*84, 95+i*84, 97+i*84, 81+i*84])
                    outlet_faces.append([95+i*84, 109+i*84, 111+i*84, 97+i*84])
                solution_to_air_faces = [[10, 24, 25, 11]]
                for i in range(geometry.n_wires):
                    solution_to_air_faces.append([
                        24+i*84, 38+i*84, 39+i*84, 25+i*84])
                    solution_to_air_faces.append([
                        38+i*84, 60+i*84, 61+i*84, 39+i*84])
                    solution_to_air_faces.append([
                        60+i*84, 78+i*84, 79+i*84, 61+i*84])
                    solution_to_air_faces.append([
                        78+i*84, 94+i*84, 95+i*84, 79+i*84])
                    solution_to_air_faces.append([
                        94+i*84, 108+i*84, 109+i*84, 95+i*84])
                wall_faces = [[12, 26, 27, 13]]
                for i in range(geometry.n_wires):
                    wall_faces.append([26+i*84, 40+i*84, 41+i*84, 27+i*84])
                    wall_faces.append([40+i*84, 62+i*84, 63+i*84, 41+i*84])
                    wall_faces.append([62+i*84, 80+i*84, 81+i*84, 63+i*84])
                    wall_faces.append([80+i*84, 96+i*84, 97+i*84, 81+i*84])
                    wall_faces.append([96+i*84, 110+i*84, 111+i*84, 97+i*84])
                front_back_faces = [
                    [10, 12, 13, 11],
                        [24+geometry.n_wires*84,
                        26+geometry.n_wires*84,
                        27+geometry.n_wires*84,
                        25+geometry.n_wires*84]]

                boundary = [
                    'inlet',
                    {'type': 'patch', 'faces': inlet_faces},
                    'outlet',
                    {'type': 'patch', 'faces': outlet_faces},
                    'wall',
                    {'type': 'patch',
                    'faces': wall_faces},
                    'front_back',
                    {'type': 'patch',
                    'faces': front_back_faces},
                    'solution_to_air',
                    {'type': 'patch',
                    'faces': solution_to_air_faces}]
            else:
                print(f'ERROR: Region {region} not defined.')
                sys.exit(1)

            geometry.edges[region] = edges
            geometry.blocks[region] = blocks
            geometry.boundary[region] = boundary
            geometry.mergePatchPairs[region] = mergePatchPairs
