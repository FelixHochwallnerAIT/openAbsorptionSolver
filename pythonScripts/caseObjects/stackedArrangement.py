#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12.05.2023

Felix Hochwallner

Class definition stacked arrangement simulations, to view, edit and simulate 
stacked arrangement OpenFOAM cases.
"""

# Import standard python modules
import os
import numpy as np
from scipy.optimize import newton
import copy

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.ParsedBlockMeshDict import ParsedBlockMeshDict
from PyFoam.Basics.DataStructures import Vector as OFVector
from PyFoam.Basics.DataStructures import Field as OFField

from Ofpp import parse_internal_field, parse_field_all, parse_boundary_field

from .baseCaseOpenFoam import baseCaseOpenFoam

class stackedArrangement(baseCaseOpenFoam):
    """
    Stores an OpenFoam stacked arrangement case.
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
        self.classtype = 'stackedArrangement'
        self.solver.multiRegion = True
        self.solver.airRegions = ['air'] 
        self.solver.solutionRegions = (
            ['solution%i' %(i+1) for i in range(self.geometry.n_steps)])
        self.solver.regions = (
            self.solver.airRegions + self.solver.solutionRegions)
        self.setFilmThickness()  
        self.initalize(case=case)
        self.createFolders()  
        self.prepareInputFiles()
        self.writeControlDict()

    def writeControlDict(self):
        """
        Writes the controlDict functions
        """  
        print("Writing the control dict functions ....")
        controlDict = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'system',
                'controlDict'))
        functions = controlDict['functions']
        # Create residuals function for every solution region
        for region in self.solver.solutionRegions:
            functions[f'residuals_{region}'] = copy.deepcopy(
                functions['residuals_solution1'])
            functions[f'residuals_{region}']['region'] = region
        # Create air_to_solution patches list
        patches = ['air_to_'+region for region in self.solver.solutionRegions]
        # Add patches to functions
        functions['absorbedMassFluxAir']['patches'] = patches
        functions['absorbedMassFluxAirObject']['patches'] = patches
        functions['sensibleHeatFluxAir']['patches'] = patches
        functions['sensibleHeatFluxDensityAirObject']['patches'] = patches
        # Write
        controlDict.writeFile()
        print("control dict functions written!")

    def createFolders(self):
        """
        Creates the necessary folders
        """    
        print("Creating region folders ....")
        # Write regionProperties file
        regionProperties = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'constant',
                'regionProperties'))
        regionProperties['regions'][1] = self.solver.regions
        regionProperties['regions'][3] = (
            ['true' for region in self.solver.airRegions] 
            + ['false' for region in self.solver.solutionRegions])
        regionProperties.writeFile()
        # Create region folders
        for region in self.solver.solutionRegions:
            # Create 0 folders if not existant
            folder_0 = os.path.join(
                self.solver.case_folder,
                '0')
            if region not in os.listdir(folder_0):
                os.system(
                    f'cp -r '
                    + os.path.join(folder_0, 'solution1')
                    + ' '
                    + os.path.join(folder_0, region))
            # Create constant folder if not existant
            folder_constant = os.path.join(
                self.solver.case_folder,
                'constant')
            if region not in os.listdir(folder_constant):
                os.system(
                    f'cp -r '
                    + os.path.join(folder_constant, 'solution1')
                    + ' '
                    + os.path.join(folder_constant, region))
            # Create system folder if not existant
            folder_system = os.path.join(
                self.solver.case_folder,
                'system')
            if region not in os.listdir(folder_system):
                os.system(
                    f'cp -r '
                    + os.path.join(folder_system, 'solution1')
                    + ' '
                    + os.path.join(folder_system, region))
        print("Region folders created!")

    def prepareInputFiles(self):
        """
        Prepares the input files cw, p, pw, T, U
        """
        print("Preparing input files ....")
        folder_0 = os.path.join(
                self.solver.case_folder,
                '0')
        # air regions
        # cw
        for region in self.solver.airRegions:
            cw = ParsedParameterFile(
                os.path.join(
                    folder_0,
                    region,
                    'cw'
                )
            )
            i = 0
            for region in self.solver.solutionRegions:
                cw['boundaryField'][f'air_to_{region}'] = {
                    'type' : 'regionCoupleVaporPressure',
                    'neighbourRegionName': region,
                    'neighbourPatchName': 'solution_to_air',
                    'neighbourFieldName': 'pw',
                    'relaxationFactor': 0.5,
                    'value': '$internalField'
                }
                cw['boundaryField']['symmetry%s' %(region[8:])] = {
                    'type': 'symmetry'
                }
                if i > 0:
                    cw['boundaryField']['symmetry_gap1_%i' %i] = {
                        'type': 'symmetry'
                    }
                    cw['boundaryField']['symmetry_gap2_%i' %i] = {
                        'type': 'symmetry'
                    }
                i += 1
            cw['boundaryField']['symmetry_inlet1'] = {
                'type': 'symmetry'
            }
            cw['boundaryField']['symmetry_inlet2'] = {
                'type': 'symmetry'
            }
            cw['boundaryField']['symmetry_outlet1'] = {
                'type': 'symmetry'
            }
            cw['boundaryField']['symmetry_outlet2'] = {
                'type': 'symmetry'
            }
            cw.writeFile()
        # p
        for region in self.solver.airRegions:
            p = ParsedParameterFile(
                os.path.join(
                    folder_0,
                    region,
                    'p'
                )
            )
            i = 0
            for region in self.solver.solutionRegions:
                p['boundaryField'][f'air_to_{region}'] = {
                    'type' : 'zeroGradient'
                }
                p['boundaryField']['symmetry%s' %(region[8:])] = {
                    'type': 'symmetry'
                }   
                if i > 0:
                    p['boundaryField']['symmetry_gap1_%i' %i] = {
                        'type': 'symmetry'
                    }   
                    p['boundaryField']['symmetry_gap2_%i' %i] = {
                        'type': 'symmetry'
                    }   
                i += 1
            p['boundaryField']['symmetry_inlet1'] = {
                'type': 'symmetry'
            }
            p['boundaryField']['symmetry_inlet2'] = {
                'type': 'symmetry'
            }
            p['boundaryField']['symmetry_outlet1'] = {
                'type': 'symmetry'
            }
            p['boundaryField']['symmetry_outlet2'] = {
                'type': 'symmetry'
            }
            p.writeFile()
        # pw
        for region in self.solver.airRegions:
            pw = ParsedParameterFile(
                os.path.join(
                    folder_0,
                    region,
                    'pw'
                )
            )
            i = 0
            for region in self.solver.solutionRegions:
                pw['boundaryField'][f'air_to_{region}'] = {
                    'type' : 'calculated',
                    'value': '$internalField'
                }
                pw['boundaryField']['symmetry%s' %(region[8:])] = {
                    'type': 'symmetry'
                }
                if i > 0:
                    pw['boundaryField']['symmetry_gap1_%i' %i] = {
                        'type': 'symmetry'
                    }  
                    pw['boundaryField']['symmetry_gap2_%i' %i] = {
                        'type': 'symmetry'
                    }  
                i += 1
            pw['boundaryField']['symmetry_inlet1'] = {
                'type': 'symmetry'
            }
            pw['boundaryField']['symmetry_inlet2'] = {
                'type': 'symmetry'
            }
            pw['boundaryField']['symmetry_outlet1'] = {
                'type': 'symmetry'
            }
            pw['boundaryField']['symmetry_outlet2'] = {
                'type': 'symmetry'
            }
            pw.writeFile()
        # T
        for region in self.solver.airRegions:
            T = ParsedParameterFile(
                os.path.join(
                    folder_0,
                    region,
                    'T'
                )
            )
            i = 0
            for region in self.solver.solutionRegions:
                T['boundaryField'][f'air_to_{region}'] = {
                    'type' : 'regionCoupleTemperature',
                    'neighbourRegionName': region,
                    'neighbourPatchName': 'solution_to_air',
                    'neighbourFieldName': 'T',
                    'relaxationFactor': 0.5,
                    'value': '$internalField'
                }
                T['boundaryField']['symmetry%s' %(region[8:])] = {
                    'type': 'symmetry'
                }
                if i > 0:
                    T['boundaryField']['symmetry_gap1_%i' %i] = {
                        'type': 'symmetry'
                    }  
                    T['boundaryField']['symmetry_gap2_%i' %i] = {
                        'type': 'symmetry'
                    }  
                i += 1
            T['boundaryField']['symmetry_inlet1'] = {
                'type': 'symmetry'
            }
            T['boundaryField']['symmetry_inlet2'] = {
                'type': 'symmetry'
            }
            T['boundaryField']['symmetry_outlet1'] = {
                'type': 'symmetry'
            }
            T['boundaryField']['symmetry_outlet2'] = {
                'type': 'symmetry'
            }
            T.writeFile()
        # U
        for region in self.solver.airRegions:
            U = ParsedParameterFile(
                os.path.join(
                    folder_0,
                    region,
                    'U'
                )
            )
            i = 0
            for region in self.solver.solutionRegions:
                U['boundaryField'][f'air_to_{region}'] = {
                    'type' : 'noSlip'
                }
                U['boundaryField']['symmetry%s' %(region[8:])] = {
                    'type': 'symmetry'
                }
                if i > 0:
                    U['boundaryField']['symmetry_gap1_%i' %i] = {
                        'type': 'symmetry'
                    }  
                    U['boundaryField']['symmetry_gap2_%i' %i] = {
                        'type': 'symmetry'
                    }  
                i += 1
            U['boundaryField']['symmetry_inlet1'] = {
                'type': 'symmetry'
            }
            U['boundaryField']['symmetry_inlet2'] = {
                'type': 'symmetry'
            }
            U['boundaryField']['symmetry_outlet1'] = {
                'type': 'symmetry'
            }
            U['boundaryField']['symmetry_outlet2'] = {
                'type': 'symmetry'
            }
            U.writeFile()
        # solution regions
        # cw
        for region in self.solver.solutionRegions:
            cw = ParsedParameterFile(
                os.path.join(
                    folder_0,
                    region,
                    'cw'
                )
            )
            cw['boundaryField']['solution_to_air']['neighbourPatchName'] = (
                f'air_to_{region}'
            )
            cw.writeFile()
        # T
        for region in self.solver.solutionRegions:
            T = ParsedParameterFile(
                os.path.join(
                    folder_0,
                    region,
                    'T'
                )
            )
            T['boundaryField']['solution_to_air']['neighbourPatchName'] = (
                f'air_to_{region}'
            )
            T.writeFile()
        print("Input files prepared!")

    def calculateVertices(self):
        """
        Calculates the vertices of the blockMesh
        """
        geometry = self.geometry
        # Calculate stuff
        geometry.x = [0]
        for i in range(geometry.n_steps):
            geometry.x += [geometry.x[-1] + geometry.l]
            if geometry.n_steps > 1:
                geometry.x += [geometry.x[-1] + geometry.x_gap]

        geometry.y = [0, 
            geometry.d_wall/2, 
            geometry.d_wall/2 + geometry.d_Nu, 
            geometry.b/2 - geometry.d_Nu - geometry.d_wall/2, 
            geometry.b/2 - geometry.d_wall/2, 
            geometry.b/2]

        geometry.z = [0, geometry.h]

        geometry.dya = (geometry.b/2) / geometry.ny_air
        geometry.nya1 = int(round((geometry.d_wall/2) / geometry.dya))
        geometry.nya2 = int(round((geometry.d_Nu) / geometry.dya))
        geometry.nya3 = int(round(
            geometry.ny_air - 2*geometry.nya1 - 2*geometry.nya2))
        assert geometry.nya2 > 10 or geometry.overwrite_check, (
            f"Less than 8 cells in solution gap, increase ny_air")
        assert geometry.nya3 > 10 or geometry.overwrite_check, (
            f"Less than 8 cells in air gap, increase ny_air")
        geometry.nx_l = int(round(geometry.l / geometry.dya))
        geometry.nx_x_gap = int(round(geometry.x_gap / geometry.dya))
        
        dx = geometry.dya
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
                print("WARNING: Could not find grading for inlet or outlet.")
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
        if geometry.calc_ny_sol:
            geometry.ny_sol = geometry.nya2

        # Number of cells
        geometry.n_cells = {'total': 0}
        for region in self.solver.regions:
            if region == 'air':
                # geometry.n_cells[region] = (geometry.n_steps 
                #     * geometry.nx * geometry.ny_air * geometry.nz)
                # Block 1
                geometry.n_cells[region] = (
                    geometry.n_steps 
                    * geometry.nx_l * geometry.nya1 * geometry.nz)
                # Block 2
                geometry.n_cells[region] += (
                    geometry.n_steps 
                    * geometry.nx_l * geometry.nya2 * geometry.nz)
                # Block 3
                geometry.n_cells[region] += (
                    geometry.n_steps 
                    * geometry.nx_l * geometry.nya3 * geometry.nz)
                # Block 4
                geometry.n_cells[region] += (
                    (geometry.n_steps-1)
                    * geometry.nx_x_gap * geometry.nya1 * geometry.nz)
                # Block 5
                geometry.n_cells[region] += (
                    (geometry.n_steps-1)
                    * geometry.nx_x_gap * geometry.nya2 * geometry.nz)
                # Block 6
                geometry.n_cells[region] += (
                    (geometry.n_steps-1)
                    * geometry.nx_x_gap * geometry.nya3 * geometry.nz)
                # Block 7
                geometry.n_cells[region] += (
                    (geometry.n_steps-1)
                    * geometry.nx_x_gap * geometry.nya2 * geometry.nz)
                # Block 8
                geometry.n_cells[region] += (
                    (geometry.n_steps-1)
                    * geometry.nx_x_gap * geometry.nya1 * geometry.nz)
                # Inlet block
                geometry.n_cells[region] += (
                    geometry.nx_inlet 
                    * (2*geometry.nya1 + 2*geometry.nya2 + geometry.nya3)
                    * geometry.nz)
                # Outlet block
                geometry.n_cells[region] += (
                    geometry.nx_outlet
                    * (2*geometry.nya1 + 2*geometry.nya2 + geometry.nya3)
                    * geometry.nz)
            elif 'solution' in region:
                geometry.n_cells[region] = (
                    geometry.nx_l * geometry.ny_sol * geometry.nz)
            else:
                raise ValueError(f"Region {region} not recognized.")
            geometry.n_cells['total'] += geometry.n_cells[region]

        # Calculate vertices
        geometry.vertices = []
        for i in range(geometry.n_steps*2):
            geometry.vertices += [
            OFVector(geometry.x[i], geometry.y[5], geometry.z[0]),
            OFVector(geometry.x[i], geometry.y[5], geometry.z[1]),
            OFVector(geometry.x[i], geometry.y[4], geometry.z[0]),
            OFVector(geometry.x[i], geometry.y[4], geometry.z[1]),
            OFVector(geometry.x[i], geometry.y[3], geometry.z[0]),
            OFVector(geometry.x[i], geometry.y[3], geometry.z[1]),
            OFVector(geometry.x[i], geometry.y[2], geometry.z[0]),
            OFVector(geometry.x[i], geometry.y[2], geometry.z[1]),
            OFVector(geometry.x[i], geometry.y[1], geometry.z[0]),
            OFVector(geometry.x[i], geometry.y[1], geometry.z[1]),
            OFVector(geometry.x[i], geometry.y[0], geometry.z[0]),
            OFVector(geometry.x[i], geometry.y[0], geometry.z[1])
            ]
        vertices = len(geometry.vertices)
        # Inlet block
        for i in range(12):
            geometry.vertices += [
                OFVector(
                    geometry.vertices[i][0] - geometry.x_inlet,
                    geometry.vertices[i][1],
                    geometry.vertices[i][2])]
        # Outlet block
        for i in range(12):
            geometry.vertices += [
                OFVector(
                    geometry.vertices[i+vertices-12][0] + geometry.x_outlet,
                    geometry.vertices[i+vertices-12][1],
                    geometry.vertices[i+vertices-12][2])]

        # Create emptpy blockMeshDict dicts
        geometry.edges = {}
        geometry.blocks = {}
        geometry.boundary = {}
        geometry.mergePatchPairs = {}

    def calculateAirBlockMeshDict(self):
        """
        Calculates the blockMeshDict of the air regions
        """
        # Create edges
        edges = []
        # Create blocks
        blocks = []
        for i in range(self.geometry.n_steps):
            shift = 48 * int(i/2)
            block1 = [10, 22, 20, 8, 11, 23, 21, 9]
            block2 = [8, 20, 18, 6, 9, 21, 19, 7]
            block3 = [6, 18, 16, 4, 7, 19, 17, 5]
            # gap
            block4 = [22, 34, 32, 20, 23, 35, 33, 21]
            block5 = [20, 32, 30, 18, 21, 33, 31, 19]
            block6 = [18, 30, 28, 16, 19, 31, 29, 17]
            block7 = [16, 28, 26, 14, 17, 29, 27, 15]
            block8 = [14, 26, 24, 12, 15, 27, 25, 13]
            # block 1
            blocks.append('hex')
            if not i%2: # even
                blocks.append([x+shift for x in block1])
                blocks.append(
                    OFVector(
                        self.geometry.nx_l, 
                        self.geometry.nya1, 
                        self.geometry.nz))
            else:
                blocks.append([x+20+shift for x in block1])
                blocks.append(
                    OFVector(
                        self.geometry.nx_l, 
                        self.geometry.nya3, 
                        self.geometry.nz))
            blocks.append('simpleGrading')
            blocks.append(OFVector(1, 1, 1))
            # block 2
            blocks.append('hex')
            if not i%2: # even
                blocks.append([x+shift for x in block2])
            else:
                blocks.append([x+20+shift for x in block2])
            blocks.append(
                OFVector(
                    self.geometry.nx_l, 
                    self.geometry.nya2, 
                    self.geometry.nz))
            blocks.append('simpleGrading')
            blocks.append(OFVector(1, 1, 1))
            # block 3
            blocks.append('hex')
            if not i%2: # even
                blocks.append([x+shift for x in block3])
                blocks.append(
                    OFVector(
                        self.geometry.nx_l, 
                        self.geometry.nya3, 
                        self.geometry.nz))
            else:
                blocks.append([x+20+shift for x in block3])
                blocks.append(
                    OFVector(
                        self.geometry.nx_l, 
                        self.geometry.nya1, 
                        self.geometry.nz))
            blocks.append('simpleGrading')
            blocks.append(OFVector(1, 1, 1))
            if i > 0:
                shift_gap = int(24 * (i-1))
                # block 4
                blocks.append('hex')
                blocks.append([x+shift_gap for x in block4])
                blocks.append(
                    OFVector(
                        self.geometry.nx_x_gap,
                        self.geometry.nya1,
                        self.geometry.nz))
                blocks.append('simpleGrading')
                blocks.append(OFVector(1, 1, 1))
                # block 5
                blocks.append('hex')
                blocks.append([x+shift_gap for x in block5])
                blocks.append(
                    OFVector(
                        self.geometry.nx_x_gap,
                        self.geometry.nya2,
                        self.geometry.nz))
                blocks.append('simpleGrading')
                blocks.append(OFVector(1, 1, 1))
                # block 6
                blocks.append('hex')
                blocks.append([x+shift_gap for x in block6])
                blocks.append(
                    OFVector(
                        self.geometry.nx_x_gap,
                        self.geometry.nya3,
                        self.geometry.nz))
                blocks.append('simpleGrading')
                blocks.append(OFVector(1, 1, 1))
                # block 7
                blocks.append('hex')
                blocks.append([x+shift_gap for x in block7])
                blocks.append(
                    OFVector(        
                        self.geometry.nx_x_gap,
                        self.geometry.nya2,
                        self.geometry.nz))
                blocks.append('simpleGrading')
                blocks.append(OFVector(1, 1, 1))
                # block 8
                blocks.append('hex')
                blocks.append([x+shift_gap for x in block8])
                blocks.append(
                    OFVector(
                        self.geometry.nx_x_gap,
                        self.geometry.nya1,
                        self.geometry.nz))
                blocks.append('simpleGrading')
                blocks.append(OFVector(1, 1, 1))
        # inlet block
        vertices = len(self.geometry.vertices)-24
        blocks.append('hex')
        blocks.append([
            vertices + 10,
            10,
            8,
            vertices + 8,
            vertices + 11,
            11,
            9,
            vertices + 9])
        blocks.append(
            OFVector(
                self.geometry.nx_inlet,
                self.geometry.nya1,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_inlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        blocks.append('hex')
        blocks.append([
            vertices + 8,
            8,
            6,
            vertices + 6,
            vertices + 9,
            9,
            7,
            vertices + 7])
        blocks.append(
            OFVector(
                self.geometry.nx_inlet,
                self.geometry.nya2,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_inlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        blocks.append('hex')
        blocks.append([
            vertices + 6,
            6,
            4,
            vertices + 4,
            vertices + 7,
            7,
            5,
            vertices + 5])
        blocks.append(
            OFVector(
                self.geometry.nx_inlet,
                self.geometry.nya3,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_inlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        blocks.append('hex')
        blocks.append([
            vertices + 4,
            4,
            2,
            vertices + 2,
            vertices + 5,
            5,
            3,
            vertices + 3])
        blocks.append(
            OFVector(
                self.geometry.nx_inlet,
                self.geometry.nya2,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_inlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        blocks.append('hex')
        blocks.append([
            vertices + 2,
            2,
            0,
            vertices + 0,
            vertices + 3,
            3,
            1,
            vertices + 1])
        blocks.append(
            OFVector(
                self.geometry.nx_inlet,
                self.geometry.nya1,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_inlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        # outlet block
        blocks.append('hex')
        blocks.append([
            vertices - 2,
            vertices + 22,
            vertices + 20,
            vertices - 4,
            vertices - 1,
            vertices + 23,
            vertices + 21,
            vertices - 3])
        blocks.append(
            OFVector(
                self.geometry.nx_outlet,
                self.geometry.nya1,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_outlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        blocks.append('hex')
        blocks.append([
            vertices - 4,
            vertices + 20,
            vertices + 18,
            vertices - 6,
            vertices - 3,
            vertices + 21,
            vertices + 19,
            vertices - 5])
        blocks.append(
            OFVector(
                self.geometry.nx_outlet,
                self.geometry.nya2,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_outlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        blocks.append('hex')
        blocks.append([
            vertices - 6,
            vertices + 18,
            vertices + 16,
            vertices - 8,
            vertices - 5,
            vertices + 19,
            vertices + 17,
            vertices - 7])
        blocks.append(
            OFVector(
                self.geometry.nx_outlet,
                self.geometry.nya3,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_outlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        blocks.append('hex')
        blocks.append([
            vertices - 8,
            vertices + 16,
            vertices + 14,
            vertices - 10,
            vertices - 7,
            vertices + 17,
            vertices + 15,
            vertices - 9])
        blocks.append(
            OFVector(
                self.geometry.nx_outlet,
                self.geometry.nya2,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_outlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        blocks.append('hex')
        blocks.append([
            vertices - 10,
            vertices + 14,
            vertices + 12,
            vertices - 12,
            vertices - 9,
            vertices + 15,
            vertices + 13,
            vertices - 11])
        blocks.append(
            OFVector(
                self.geometry.nx_outlet,
                self.geometry.nya1,
                self.geometry.nz))
        blocks.append('simpleGrading')
        if self.geometry.grading_inlet_outlet:
            blocks.append(OFVector(self.geometry.grading_outlet, 1, 1))
        else:
            blocks.append(OFVector(1, 1, 1))
        
        # Create boundaries
        boundary = []
        # inlet
        faces_inlet = [
            [10, 8, 9, 11], 
            [8, 6, 7, 9], 
            [6, 4, 5, 7],
            [4, 2, 3, 5],
            [2, 0, 1, 3]]
        faces_inlet_new = []
        for i in range(5):
            faces_inlet_new.append(
                [x+vertices for x in faces_inlet[i]])
        boundary.append('inlet')
        boundary.append(
            {'type': 'patch',
            'faces': faces_inlet_new})
        # outlet
        faces_outlet = []
        for face in faces_inlet:
            if self.geometry.n_steps > 1:
                faces_outlet.append(
                    [x+24*self.geometry.n_steps-12+24 for x in face])
            else:
                faces_outlet.append(
                    [x+36 for x in face])
        boundary.append('outlet')
        boundary.append(
            {'type': 'patch',
            'faces': faces_outlet})
        # top_bottom
        faces_top_bottom_first_step = [
            [10, 22, 20, 8],
            [8, 20, 18, 6],
            [6, 18, 16, 4],
            [11, 23, 21, 9],
            [9, 21, 19, 7],
            [7, 19, 17, 5]]
        faces_top_bottom = []
        for i in range(self.geometry.n_steps):
            shift = 24*i
            for face in faces_top_bottom_first_step:
                if not i%2: # even step
                    faces_top_bottom.append(
                        [x+shift for x in face])
                else: # uneven step
                    faces_top_bottom.append(
                        [x+shift-4 for x in face])
        for i in range(1, self.geometry.n_steps):
            shift = 24*(i-1)
            faces_top_bottom.append(
                [22+shift, 34+shift, 32+shift, 20+shift])
            faces_top_bottom.append(
                [20+shift, 32+shift, 30+shift, 18+shift])
            faces_top_bottom.append(
                [18+shift, 30+shift, 28+shift, 16+shift])
            faces_top_bottom.append(
                [16+shift, 28+shift, 26+shift, 14+shift])
            faces_top_bottom.append(
                [14+shift, 26+shift, 24+shift, 12+shift])
            faces_top_bottom.append(
                [23+shift, 35+shift, 33+shift, 21+shift])
            faces_top_bottom.append(
                [21+shift, 33+shift, 31+shift, 19+shift])
            faces_top_bottom.append(
                [19+shift, 31+shift, 29+shift, 17+shift])
            faces_top_bottom.append(
                [17+shift, 29+shift, 27+shift, 15+shift])
            faces_top_bottom.append(
                [15+shift, 27+shift, 25+shift, 13+shift])
        # inlet block
        faces_top_bottom.append([
            vertices + 10,
            10,
            8,
            vertices + 8])
        faces_top_bottom.append([
            vertices + 8,
            8,
            6,
            vertices + 6])
        faces_top_bottom.append([
            vertices + 6,
            6,
            4,
            vertices + 4])
        faces_top_bottom.append([
            vertices + 4,
            4,
            2,
            vertices + 2])
        faces_top_bottom.append([
            vertices + 2,
            2,
            0,
            vertices + 0])
        faces_top_bottom.append([
            vertices + 11,
            11,
            9,
            vertices + 9])
        faces_top_bottom.append([
            vertices + 9,
            9,
            7,
            vertices + 7])
        faces_top_bottom.append([
            vertices + 7,
            7,
            5,
            vertices + 5])
        faces_top_bottom.append([
            vertices + 5,
            5,
            3,
            vertices + 3])
        faces_top_bottom.append([
            vertices + 3,
            3,
            1,
            vertices + 1])
        # outlet block
        faces_top_bottom.append([
            vertices - 2,
            vertices + 22,
            vertices + 20,
            vertices - 4])
        faces_top_bottom.append([
            vertices - 4,
            vertices + 20,
            vertices + 18,
            vertices - 6])
        faces_top_bottom.append([
            vertices - 6,
            vertices + 18,
            vertices + 16,
            vertices - 8])
        faces_top_bottom.append([
            vertices - 8,
            vertices + 16,
            vertices + 14,
            vertices - 10])
        faces_top_bottom.append([
            vertices - 10,
            vertices + 14,
            vertices + 12,
            vertices - 12])
        faces_top_bottom.append([
            vertices - 1,
            vertices + 23,
            vertices + 21,
            vertices - 3])
        faces_top_bottom.append([
            vertices - 3,
            vertices + 21,
            vertices + 19,
            vertices - 5])
        faces_top_bottom.append([
            vertices - 5,
            vertices + 19,
            vertices + 17,
            vertices - 7])
        faces_top_bottom.append([
            vertices - 7,
            vertices + 17,
            vertices + 15,
            vertices - 9])
        faces_top_bottom.append([
            vertices - 9,
            vertices + 15,
            vertices + 13,
            vertices - 11])
        boundary.append('top_bottom')
        boundary.append(
            {'type': 'patch',
            'faces': faces_top_bottom})
        # air_to_solution
        faces_air_to_solution_first_step = [4, 16, 17, 5]
        for i in range(1,self.geometry.n_steps+1):
            if not i%2: # even step
                shift = 24*(i-1)+2
            else: # uneven step
                shift = 24*(i-1)
            faces_air_to_solution = [[
                x+shift for x in faces_air_to_solution_first_step]]
            boundary.append('air_to_solution%i' %i)
            boundary.append(
                {'type': 'patch',
                'faces': faces_air_to_solution})
        # symmetry
        faces_symmetry_first_step = [10, 22, 23, 11]
        for i in range(1,self.geometry.n_steps+1):
            if not i%2: # even step
                shift = 24*(i-1)-10
            else: # uneven step
                shift = 24*(i-1)
            faces_symmetry = [[
                x+shift for x in faces_symmetry_first_step]]
            boundary.append('symmetry%i' %i)
            boundary.append(
                {'type': 'symmetry',
                'faces': faces_symmetry})
        # symmetry gap 1
        faces_symmetry_first_gap1 = [22, 34, 35, 23]
        for i in range(1,self.geometry.n_steps):
            faces_symmetry_gap1 = [
                [x+24*(i-1) for x in faces_symmetry_first_gap1]]
            boundary.append('symmetry_gap1_%i' %i)
            boundary.append(
                {'type': 'symmetry',
                'faces': faces_symmetry_gap1})
        # symmetry gap 2
        faces_symmetry_first_gap2 = [12, 24, 25, 13]
        for i in range(1,self.geometry.n_steps):
            faces_symmetry_gap2 = [
                [x+24*(i-1) for x in faces_symmetry_first_gap2]]
            boundary.append('symmetry_gap2_%i' %i)
            boundary.append(
                {'type': 'symmetry',
                'faces': faces_symmetry_gap2})
        # symmetry inlet 1
        faces_symmetry_first_inlet1 = [
            vertices + 10,
            10,
            11,
            vertices + 11]
        boundary.append('symmetry_inlet1')
        boundary.append(
            {'type': 'symmetry',
            'faces': [faces_symmetry_first_inlet1]})
        # symmetry inlet 2
        faces_symmetry_first_inlet2 = [
            vertices + 0,
            0,
            1,
            vertices + 1]
        boundary.append('symmetry_inlet2')
        boundary.append(
            {'type': 'symmetry',
            'faces': [faces_symmetry_first_inlet2]})
        # symmetry outlet 1
        faces_symmetry_first_outlet1 = [
            vertices - 2,
            vertices + 22,
            vertices + 23,
            vertices - 1]
        boundary.append('symmetry_outlet1')
        boundary.append(
            {'type': 'symmetry',
            'faces': [faces_symmetry_first_outlet1]})
        # symmetry outlet 2
        faces_symmetry_first_outlet2 = [
            vertices - 12,
            vertices + 12,
            vertices + 13,
            vertices - 11]
        boundary.append('symmetry_outlet2')
        boundary.append(
            {'type': 'symmetry',
            'faces': [faces_symmetry_first_outlet2]})
        # wall
        faces_wall = []
        faces_wall1 = [34, 32, 33, 35]
        faces_wall2 = [32, 30, 31, 33]
        faces_wall3 = [16, 14, 15, 17]
        faces_wall4 = [14, 12, 13, 15]
        faces_wall_inlet = [0, 2, 3, 1]
        faces_solution_inlet = [2, 4, 5, 3]
        if not self.geometry.n_steps%2: # even amount of steps
            faces_wall_outlet = [
                vertices-4, vertices-2, vertices-1, vertices-3]
            faces_solution_outlet = [
                vertices-6, vertices-4, vertices-3, vertices-5]
        else: # uneven amount of steps
            faces_wall_outlet = [
                vertices-12, vertices-10, vertices-9, vertices-11]
            faces_solution_outlet = [
                vertices-10, vertices-8, vertices-7, vertices-9]
        for i in range(2,self.geometry.n_steps+1):
            shift = 24*(i-2)
            if not i%2: # even step
                faces_wall.append([x+shift for x in faces_wall1])
                faces_wall.append([x+shift for x in faces_wall2])
                faces_wall.append([x+shift for x in faces_wall3])
                faces_wall.append([x+shift for x in faces_wall4])
            else: # uneven step
                faces_wall.append([x+shift-6 for x in faces_wall1])
                faces_wall.append([x+shift-6 for x in faces_wall2])
                faces_wall.append([x+shift+6 for x in faces_wall3])
                faces_wall.append([x+shift+6 for x in faces_wall4])
        faces_wall.append(faces_wall_inlet)
        faces_wall.append(faces_solution_inlet)
        faces_wall.append(faces_wall_outlet)
        faces_wall.append(faces_solution_outlet)
        boundary.append('wall')
        boundary.append(
            {'type': 'patch',
            'faces': faces_wall})

        # Create mergePatchPairs
        mergePatchPairs = []

        for region in self.solver.airRegions:
            self.geometry.edges[region] = edges
            self.geometry.blocks[region] = blocks
            self.geometry.boundary[region] = boundary
            self.geometry.mergePatchPairs[region] = mergePatchPairs

    def calculateSolutionBlockMeshDict(self):
        """
        Calculates the blockMeshDict of the air region
        """
        # Create edges
        edges = []

        # Create mergePatchPairs
        mergePatchPairs = []

        for region in self.solver.solutionRegions:
            # Which solution region?
            i = int(region[8:])
            if not i%2: # even
                shift = 24*(i-1)+4
            else: # uneven
                shift = 24*(i-1)

            # Create blocks
            blocks = []
            block1 = [4, 16, 14, 2, 5, 17, 15, 3]
            # block 1
            blocks.append('hex')
            blocks.append([x+shift for x in block1])
            blocks.append(
                OFVector(
                    self.geometry.nx_l, 
                    self.geometry.ny_sol, 
                    self.geometry.nz))
            blocks.append('simpleGrading')
            if not i%2: # even
                blocks.append(
                    OFVector(1, 1/self.geometry.solution.grading, 1))
            else: # uneven
                blocks.append(
                    OFVector(1, self.geometry.solution.grading, 1))

            # Create boundaries
            boundary = []
            # inlet
            faces_inlet_first = [4, 16, 14, 2]
            faces_inlet = [[x+shift for x in faces_inlet_first]]
            boundary.append('inlet')
            boundary.append(
                {'type': 'patch',
                'faces': faces_inlet})
            # outlet
            faces_outlet_first = [5, 17, 15, 3]
            faces_outlet = [[x+shift for x in faces_outlet_first]]
            boundary.append('outlet')
            boundary.append(
                {'type': 'patch',
                'faces': faces_outlet})
            # wall
            faces_wall_first = [2, 14, 15, 3]
            if not i%2: # even
                faces_wall = [[x+shift+2 for x in faces_wall_first]]
            else: # uneven
                faces_wall = [[x+shift for x in faces_wall_first]]
            boundary.append('wall')
            boundary.append(
                {'type': 'patch',
                'faces': faces_wall})
            # front_back
            faces_front_back_first = [
                [4, 2, 3, 5],
                [16, 14, 15, 17]]
            faces_front_back = []
            for face in faces_front_back_first:
                faces_front_back.append([x+shift for x in face])
            boundary.append('front_back')
            boundary.append(
                {'type': 'patch',
                'faces': faces_front_back})
            # solution_to_air
            faces_solution_to_air_first = [4, 16, 17, 5]
            if not i%2: # even
                faces_solution_to_air = [
                    [x+shift-2 for x in faces_solution_to_air_first]]
            else: # uneven
                faces_solution_to_air = [
                    [x+shift for x in faces_solution_to_air_first]]
            boundary.append('solution_to_air')
            boundary.append(
                {'type': 'patch',
                'faces': faces_solution_to_air})

            self.geometry.edges[region] = edges
            self.geometry.blocks[region] = blocks
            self.geometry.boundary[region] = boundary
            self.geometry.mergePatchPairs[region] = mergePatchPairs

    def set2D(self):
        """
        Sets the case 2D
        """
        # Set region properties
        regionProperties = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'constant',
                'regionProperties'))
        regionProperties['regions'][1] = ['air']
        regionProperties['regions'][3] = ['true']
        regionProperties.writeFile()
        # Set blockMeshDict top_bottom to empty
        blockMeshDict = ParsedBlockMeshDict(
            os.path.join(
                self.solver.case_folder,
                'system',
                'air',
                'blockMeshDict'))
        idx_top_bot = blockMeshDict['boundary'].index('top_bottom')
        blockMeshDict['boundary'][idx_top_bot+1]['type'] = 'empty'
        blockMeshDict.writeFile()
        # Create new mesh
        self.createMesh()
        # Set boundaries to empty
        fields = ['U', 'p', 'T', 'cw', 'pw']
        for field in fields:
            field_air = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    '0',
                    'air',
                    field))
            field_air['boundaryField']['top_bottom']['type'] = 'empty'
            field_air.writeFile()
        # Set cw field
        # Calculate C_a_int
        if self.material.solution.medium == 'LiBr':
            from ..fluidProperties import LiBrProperties \
                as desiccantProperties
        else:
            raise ValueError(
                'Medium not implemented!')
        p_s_int = desiccantProperties.pw(
            1 - self.initialization.solution.C_in,
            self.initialization.solution.T_wall)
        C_a_int = humidAirProperties.cwPw(
            1e5,
            p_s_int)
        cw = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    '0',
                    'air',
                    'cw'))
        for region in self.solver.solutionRegions:
            cw['boundaryField'][f'air_to_{region}'] = {
                'type': 'fixedValue',
                'value': f'uniform {C_a_int}'
            }
        cw.writeFile()
        # Set T field
        T = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    '0',
                    'air',
                    'T'))
        for region in self.solver.solutionRegions:
            T['boundaryField'][f'air_to_{region}'] = {
                'type': 'fixedValue',
                'value': ('uniform ' 
                    + str(self.initialization.solution.T_wall+273.15))
            }
        T.writeFile()
        # Delete controlDict functions
        controlDict = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'system',
                'controlDict'))
        functions = controlDict['functions']
        # Create air_to_solution patches list
        residuals = [
            'residuals_'+region for region in self.solver.solutionRegions]
        for residual in residuals:
            del functions[residual]
        controlDict.writeFile()
