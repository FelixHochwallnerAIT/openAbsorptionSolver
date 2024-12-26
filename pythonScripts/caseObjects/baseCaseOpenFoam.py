#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:17:21 2020

@author: Felix Hochwallner

Class definition baseCaseOpenFoam, for OpenFoam cases.
"""

import os
import numpy as np
import glob
import pickle
from datetime import datetime
import json

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.ParsedBlockMeshDict import ParsedBlockMeshDict
from PyFoam.Basics.DataStructures import DictProxy, TupleProxy
from PyFoam.Basics.DataStructures import Vector as OFVector
from PyFoam.Basics.DataStructures import Field as OFField

from Ofpp import parse_internal_field, parse_field_all, parse_boundary_field

from .fluidProperties import humidAirProperties

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

class baseCaseOpenFoam:
    """
    Stores an OpenFoam case.
    """
    def __init__(self, solver, geometry, material, initialization):
        """
        Constructor class
        Input:
        solver: Solver settings (solverSettings Class)
        geometry: Geometry settings (geometrySettings Class)
        material: Material settings (materialSettings Class)
        initialization: Initialization (initializationSettings Class)
        """
        self.solver = solver
        self.geometry = geometry
        self.material = material
        self.initialization = initialization
        
        self.results = dict()
        self.classtype = 'baseCase'
        self.setWriteFormat('ascii')

    def initalize(self, case=''):
        """
        Initalizes the case
        """
        print("Initalizing case ....")
        now = datetime.now()
        current_time = now.strftime("%Y-%m-%d_%H:%M:%S")
        if case != '':
            workingFolderName = case
        else:
            workingFolderName = (self.classtype+'_'+current_time)
        self.solver.working_folder_name = workingFolderName
        self.solver.working_folder = os.path.abspath(
            os.path.join(
                self.solver.working_folder,
                case))
        # Create working folder
        os.system('mkdir %s' %self.solver.working_folder)
        # Copy case folder to working folder
        os.system('cp -r %s %s' 
            %(  
                os.path.join(
                    self.solver.case_folder, 
                    '0'), 
                self.solver.working_folder))
        os.system('cp -r %s %s' 
            %(  
                os.path.join(
                    self.solver.case_folder, 
                    'constant'), 
                self.solver.working_folder))
        os.system('cp -r %s %s' 
            %(  
                os.path.join(
                    self.solver.case_folder, 
                    'system'), 
                self.solver.working_folder))
        os.system('touch %s' 
            %os.path.join(
                self.solver.working_folder, 
                '%s.foam' %workingFolderName))
        # Set case folder to working folder
        self.solver.case_folder_original = self.solver.case_folder
        self.solver.case_folder = self.solver.working_folder
        # Write json file of object parameters
        self.writeParametersJSON()
        # Create pickle of object
        pickle.dump(
            self, 
            open(
                os.path.join(
                    self.solver.case_folder,
                    "pythonCase.p"), 
                "wb"))
        print('Case initalized!')

    def setInitialization(self):
        """
        Sets the initialization settings as defined in self.initialization
        """
        print("Setting initialization settings for air region ....")
        # Set initialization for U in air region
        U_air = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                str(self.solver.startTime),
                'air',
                'U'))
        U_air['internalField'] = ('uniform ('
            + str(self.initialization.air.v_in)
            + ' 0 0)')
        U_air['boundaryField']['inlet']['type'] = 'fixedValue'
        U_air['boundaryField']['inlet']['value'] = ('uniform ('
            + str(self.initialization.air.v_in)
            + ' 0 0)')
        U_air.writeFile()
        os.system('rm '
            + os.path.join(
                self.solver.case_folder,
                str(self.solver.startTime),
                'air',
                'phi')
            + ' > /dev/null 2>&1')
        # Set initialization for T in air region
        T_air = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                str(self.solver.startTime),
                'air',
                'T'))
        T_air['internalField'] = ('uniform '
            + str(self.initialization.air.T_in+273.15))
        T_air['boundaryField']['inlet']['value'] = ('uniform '
            + str(self.initialization.air.T_in+273.15))
        T_air['boundaryField']['outlet']['inletValue'] = ('uniform '
            + str(self.initialization.air.T_in+273.15))
        for sol_reg in self.solver.solutionRegions:
            T_air['boundaryField'][f'air_to_{sol_reg}']['value'] = ('uniform '
                + str(self.initialization.solution.T_wall+273.15+5.0))
        T_air.writeFile()
        # Set initialization for cw in air region
        C_air = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                str(self.solver.startTime),
                'air',
                'cw'))
        C_air['internalField'] = ('uniform '
            + str(self.initialization.air.C_in))
        C_air['boundaryField']['inlet']['value'] = ('uniform '
            + str(self.initialization.air.C_in))
        C_air['boundaryField']['outlet']['inletValue'] = ('uniform '
            + str(self.initialization.air.C_in))
        for sol_reg in self.solver.solutionRegions:
            (C_air['boundaryField']
                [f'air_to_{sol_reg}']
                ['type']) = (
                    'regionCoupleVaporPressure')
            (C_air['boundaryField']
                [f'air_to_{sol_reg}']
                ['value']) = (
                '   $internalField')
            (C_air['boundaryField']
                [f'air_to_{sol_reg}']
                ['neighbourRegionName']) = (
                    sol_reg)
            (C_air['boundaryField']
                [f'air_to_{sol_reg}']
                ['neighbourPatchName']) = (
                    'solution_to_air')
            (C_air['boundaryField']
                [f'air_to_{sol_reg}']
                ['neighbourFieldName']) = (
                    'pw')
            (C_air['boundaryField']
                [f'air_to_{sol_reg}']
                ['relaxationFactor']) = (
                    0.5)
        C_air.writeFile()
        # Set initialization for pw in air region
        pw_air = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                str(self.solver.startTime), 
                'air', 
                'pw'))
        pw_air_inlet = humidAirProperties.pwC(
            1e5, 
            self.initialization.air.T_in, 
            self.initialization.air.C_in)
        self.initialization.air.pw_in = pw_air_inlet
        pw_air['internalField'] = f'uniform {pw_air_inlet}'
        pw_air.writeFile()
        print("Initialization of air region set!")
        # Set prescribed solution flow
        self.setPrescribedSolutionFlow()
        for region in self.solver.solutionRegions:
            print("Setting initialization settings for %s region ...." %region)
            # Set initialization for T in solution region
            T_sol = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    str(self.solver.startTime),
                    region,
                    'T'))
            T_sol['internalField'] = ('uniform '
                + str(self.initialization.solution.T_in+273.15))
            T_sol['boundaryField']['inlet']['value'] = ('uniform '
                + str(self.initialization.solution.T_in+273.15))
            T_sol['boundaryField']['wall']['type'] = 'fixedValue'
            T_sol['boundaryField']['wall']['value'] = ('uniform '
                + str(self.initialization.solution.T_wall+273.15))
            T_sol['boundaryField']['solution_to_air']['value'] = ('uniform '
                + str(self.initialization.solution.T_wall+273.15+5.0))
            T_sol.writeFile()
            # Set initialization for cw in solution region
            C_sol = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    str(self.solver.startTime),
                    region,
                    'cw'))
            C_sol['internalField'] = ('uniform '
                + str(self.initialization.solution.C_in))
            C_sol['boundaryField']['inlet']['value'] = ('uniform '
                + str(self.initialization.solution.C_in))
            C_sol.writeFile()
            # Set initialization for pw in solution region
            pw_sol = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    str(self.solver.startTime),
                    region,
                    'pw'))
            pw_sol['internalField'] = 'uniform '+str(pw_air_inlet/2)
            pw_sol.writeFile()
            print("Initialization set for %s region!" %region)

    def setFilmThickness(self):
        """
        Calculates and sets the film thickness
        """
        print("Setting film thickness ....")
        # Calculate Nusselt thickness
        g = 9.81
        delta_Nu = (
            (3.*self.material.solution.nu*self.initialization.solution.mDot)
            / (self.geometry.solution.l*self.material.solution.rho*g)
            )**(1./3)
        for geo in [self.geometry, self.geometry.solution, self.geometry.air]:
            geo.d_Nu = delta_Nu
        print("Film thickness: " + str(round(delta_Nu*1e6,2)) + "e-6 m")
        print("Film thickness set!")

    def setPrescribedSolutionFlow(self):
        """
        Sets the prescribed solution flow
        """
        def write_U_file(region):
            # Checking if phi exists, and if it exists delete it
            dir_cont = os.listdir(
                os.path.join(
                    self.solver.case_folder,
                    '0',
                    region))
            if 'phi' in dir_cont:
                os.system(f'rm {self.solver.case_folder}/0/{region}/phi')
            print("Reading C file ....")
            C_file_path = os.path.join(
                self.solver.case_folder,
                '0',
                region,
                'C')
            C_field = parse_field_all(C_file_path)
            y_cell = C_field[0][:,1]
            # Calculate Nusselt profile
            g = 9.81
            y_wall = np.mean(
                C_field[1][b'wall'][b'value'][:,1]
            )
            y_interface = np.mean(
                C_field[1][b'solution_to_air'][b'value'][:,1]
            )
            y = ((y_cell - y_wall)/(y_interface - y_wall)
                *self.geometry.d_Nu)
            u_Nu = (g/self.material.solution.nu 
                * (self.geometry.d_Nu - y/2)*y)
            print("Reading U file ....")
            U = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    '0',
                    region,
                    'U'))
            print("Reading done!")
            print("Preparing and setting U profile ....")
            U_profile = []
            for i in range(len(y_cell)):
                U_profile.append(OFVector(0,0,u_Nu[i]))
            U['internalField'] = OFField(U_profile)
            # Set inlet profile
            # Get y coordinates
            C_inlet = C_field[1][b'inlet'][b'value']
            y_cell_inlet = np.array(C_inlet[:,1])
            # Calculate Nusselt profile
            y = ((y_cell_inlet - y_wall)/(y_interface - y_wall)
                *self.geometry.d_Nu)
            u_Nu_inlet = (g/self.material.solution.nu 
                * (self.geometry.solution.d_Nu - y/2)*y)        
            U_inlet_profile = []
            for i in range(len(C_inlet)):
                U_inlet_profile.append(OFVector(0,0,u_Nu_inlet[i]))
            U['boundaryField']['inlet']['volumetricFlowRate'][-1] = (
                self.initialization.solution.mDot/self.material.solution.rho)
            U['boundaryField']['inlet']['value'] = OFField(U_inlet_profile)
            print("Preparation and setting done!")
            U.writeFile()
        for region in self.solver.solutionRegions:
            print(f"Setting prescribed Nusselt flow "
                f"for {region} region ....")
            write_U_file(region)
            print(f"Prescribed Nusselt flow "
                f"for {region} region set!")

    def setBlockMeshDict(self):
        """
        Sets the blockMeshDict.
        """
        print(f'Calculating vertices ....')
        self.calculateVertices()
        print(f'Calculaton of vertices done!')
        print(f'Calculating blockMeshDict of air region ....')
        self.calculateAirBlockMeshDict()
        print(f'Calculaton of blockMeshDict of air region done!')
        print(f'Calculating blockMeshDict of solution region ....')
        self.calculateSolutionBlockMeshDict()
        print(f'Calculaton of blockMeshDict of solution region done!')
        for region in self.solver.regions:
            print(f'Setting blockMeshDict of region {region} ....')
            blockMeshDict = ParsedBlockMeshDict(
                os.path.join(
                    self.solver.case_folder,
                    'system',
                    region,
                    'blockMeshDict'))
            blockMeshDict['vertices'] = (
                self.geometry.vertices)
            blockMeshDict['blocks'] = (
                self.geometry.blocks[region])
            blockMeshDict['edges'] = (
                self.geometry.edges[region])
            blockMeshDict['boundary'] = (
                self.geometry.boundary[region])
            blockMeshDict['mergePatchPairs'] = (
                self.geometry.mergePatchPairs[region])
            blockMeshDict.writeFile()
            print(f'blockMeshDict of region {region} set!')
        # Rewrite the JSON file
        self.writeParametersJSON()

    def writeParametersJSON(self):
        """
        Writes parameters to json file
        """
        json_object = self.toJSON()
        json_object_path = os.path.join(
            self.solver.working_folder,
            "parameters.json")
        with open(json_object_path, "w") as outfile:
            outfile.write(json_object)

    def simulate(self):
        """
        Simulates the case
        """
        # Save current working directory to go back afterwards
        current_directory = os.getcwd()
        # Got to case folder
        os.chdir(self.solver.case_folder)
        print("Starting to simulate ....")
        if self.solver.n_proc > 1:
            # Decompose and save log file, including error messages
            self.decomposeCase()
            # Simulate and save log file, including error messages
            code = os.system(
                'mpirun -np '
                +str(self.solver.n_proc)
                +' '
                +self.solver.solver
                +' -parallel > log.'
                +self.solver.solver
                +' 2>&1')
            # Reconstruct and save log file, including error messages
            if self.solver.multiRegion:
                for region in self.solver.regions:
                    if self.solver.transient:
                        os.system(
                            f'reconstructPar -region {region} '
                            f'> log.decomposePar_{region} 2>&1'
                        )
                    else:
                        os.system(
                            f'reconstructPar -latestTime -region {region} '
                            f'> log.decomposePar_{region} 2>&1'
                        )
            else:
                if self.solver.transient:
                    os.system(
                        f'reconstructPar '
                        f'> log.decomposePar 2>&1'
                    )
                else:
                    os.system(
                        f'reconstructPar -latestTime '
                        f'> log.decomposePar 2>&1'
                    )
            # Delete processor folder
            os.system("rm -r processor*")
        else:
            # Simulate and save log file, including error messages
            code = os.system(
                self.solver.solver
                +' > log.'
                +self.solver.solver
                +' 2>&1')
        # Check if simulation went well
        if code == 0:
            print("Simulation done!")
        else:
            print("ERROR: Simulation failed.")
            os.chdir(current_directory)
            raise SystemError(
                os.system(
                    'tail -1 '
                    +self.solver.case_folder
                    +'/log.'
                    +self.solver.solver))
        # Get and save latest_folder
        dir_list = os.listdir()
        result_list = []
        for dir_name in dir_list:
            if dir_name.isdigit():
                result_list.append(float(dir_name))
        latest_folder_float = max(result_list)
        latest_folder = ('%f' % latest_folder_float).rstrip('0').rstrip('.')
        self.solver.latest_folder = str(latest_folder)
        # Change back to original directory
        os.chdir(current_directory)

    def saveResults(self, foldername=''):
        print("Saving results ....")
        # Save current working directory to go back afterwards
        current_directory = os.getcwd()
        # Got to case folder
        os.chdir(self.solver.case_folder)
        # Create pickle of object
        pickle.dump(self, open("pythonCase.p", "wb"))
        # Get regions list (if regions exist)
        const_regions = glob.glob('constant/*/')
        if 'constant/polyMesh/' in const_regions:
            const_regions = []
        regions = []
        for reg in const_regions:
            regions.append(reg[9:-1])
        # Create results folder
        #os.system('mkdir -p results/'+foldername)
        results_path = os.path.join(
            self.solver.results_folder, 
            self.classtype, 
            foldername)
        os.system('mkdir -p %s' %results_path)
        # Create cell center coordinates and cell volumes in latest time folder
        for region in regions:
            os.system(
                'postProcess -func "writeCellCentres" -region %s ' %region
                + '-latestTime > log.writeCellCentres_%s' %region)
            os.system(
                'postProcess -func "writeCellVolumes" -region %s ' %region
                + '-latestTime > log.writeCellVolumes_%s' %region)
        if regions == []:
            os.system('postProcess -func "writeCellCentres"'
                    ' -latestTime > log.writeCellCentres')
            os.system('postProcess -func "writeCellVolumes"'
                    ' -latestTime > log.writeCellVolumes')
            
        # Copy 0 folder
        os.system('cp -r 0/  %s' %results_path) 
        if self.solver.transient == False:
            # Steady state case
            # Get latest result directory
            dir_list = os.listdir()
            result_list=[]
            for dir_name in dir_list:
                if dir_name.isdigit():
                    result_list.append(float(dir_name))
            latest_folder_float = max(result_list)
            latest_folder = ('%f' % latest_folder_float).rstrip('0').rstrip('.')
            self.solver.latest_folder = latest_folder
            # Copying results
            os.system('cp -r %s/ %s' %(latest_folder, results_path))
        elif self.solver.transient == True:
            # Transient Case
            # Get all directories
            # Get list of name of all folders
            folders = [name for name in os.listdir() if os.path.isdir(name)]
            # Filter only time steps
            time_folders = [name for name in folders if isfloat(name)]
            # Filter out 0 folder
            time_steps = [name for name in time_folders if name != '0']
            # Copying results
            for time in time_steps:
                os.system('cp -r %s/ %s' %(time, results_path))
        # Copying postProcessing results
        os.system('cp -r postProcessing/ %s ' %results_path +
            '> /dev/null 2>&1')
        # Copying transport properties in constant folder (if exist)
        if os.path.isfile('./constant/transportProperties'):
            os.system(
                'cp --parents constant/transportProperties %s' %results_path)
        # Copying thermophysical properties of regions (if exist)
        for region in regions:
            os.system('cp --parents constant/%s' %region +
                '/thermophysicalProperties %s' %results_path)
        # Copying blockMeshDict
        # For single region case
        if os.path.isfile('./system/blockMeshDict'):
            os.system('cp --parents system/blockMeshDict ' +
                '%s' %results_path)
        # For multi region case
        for region in regions:
            os.system('cp --parents system/%s' %region +
                '/blockMeshDict %s' %results_path)
        # Copying fvSolution
        # For single region case
        if os.path.isfile('./system/blockMeshDict'):
            os.system('cp --parents system/fvSolution ' +
                '%s' %results_path)
        # For multi region case
        for region in regions:
            os.system('cp --parents system/%s' %region +
                '/fvSolution %s' %results_path)
        # Copying fvSchemes
        # For single region case
        if os.path.isfile('./system/blockMeshDict'):
            os.system('cp --parents system/fvSchemes ' +
                '%s' %results_path)
        # For multi region case
        for region in regions:
            os.system('cp --parents system/%s' %region +
                '/fvSchemes %s' %results_path)
        # Copying fvOptions
        # For single region case
        if os.path.isfile('./system/fvOptions'):
            os.system('cp --parents system/fvOptions '+\
                '%s' %results_path)
        # For multi region case
        for region in regions:
            if os.path.isfile('./system/%s/fvOptions' %region):
                os.system('cp --parents system/%s' %region +
                    '/fvOptions %s' %results_path)
        # Copying controlDict
        if os.path.isfile('./system/controlDict'):
            os.system('cp --parents system/controlDict '+\
                '%s' %results_path)
        # Copying regionProperties
        if os.path.isfile('./constant/regionProperties'):
            os.system('cp --parents constant/regionProperties '+\
                '%s' %results_path)
        # Copying decomposeParDict
        # For single region case
        if os.path.isfile('./system/decomposeParDict'):
            os.system('cp --parents system/decomposeParDict '+\
                '%s' %results_path)
        # For multi region case
        for region in regions:
            if os.path.isfile('./system/%s/decomposeParDict' %region):
                os.system('cp --parents system/%s' %region +
                    '/decomposeParDict %s' %results_path)
        # Copying mapFieldsDict
        if os.path.isfile('./system/mapFieldsDict'):
            os.system('cp --parents system/mapFieldsDict '+\
                '%s' %results_path)
        # Save python pickle object
        os.system('cp pythonCase.p %s' %results_path)
        # Copying log-files
        os.system('cp log.* %s' %results_path)
        # Change back to original directory
        os.chdir(current_directory)
        print("Results saved!")

    def cleanCase(self):
        print("Cleaning up case ....")
        # Delete working folder
        os.system('rm -r %s > /dev/null 2>&1' %self.solver.working_folder)
        print("Case cleaned up!")

    def clearResults(self):
        print("Deleting results ....")
        # Save current working directory to go back afterwards
        current_directory = os.getcwd()
        # Got to case folder
        os.chdir(self.solver.case_folder)
        # Clean up case
        os.system('rm -r results/* > /dev/null 2>&1')
        # Change back to original directory
        os.chdir(current_directory)
        print("Results cleaned up!")

    def createMesh(self):
        """
        Creates the mesh according to the settings in blockMeshDict
        """
        print("Creating mesh ....")
        # Save current working directory to go back afterwards
        current_directory = os.getcwd()
        # Got to case folder
        os.chdir(self.solver.case_folder)
        if not self.solver.multiRegion:
            # Create mesh
            code = os.system('blockMesh 2>&1 > log.blockMesh')
            # Recreate coordinates
            code2 = os.system('postProcess -func writeCellCentres 2>&1 > log.writeCellCentres')
            # Check if everything went well
            if (code == 0) & (code2 == 0):
                print("Mesh created!")
            elif (code != 0) & (code2 == 0):
                print("ERROR: Mesh creation failed.")
                # Change back to original directory
                os.chdir(current_directory)
                raise SystemError(
                    os.system(
                        'tail -1 '
                        +os.path.join(
                            self.solver.case_folder,
                            'log.blockMesh')))
            elif (code == 0) & (code2 != 0):
                print("ERROR: Coordinate creation failed.")
                # Change back to original directory
                os.chdir(current_directory)
                raise SystemError(
                    os.system(
                        'tail -1 '
                        +os.path.join(
                            self.solver.case_folder,
                            'log.writeCellCentres')))
            else:
                print("ERROR: Both mesh and coordinate creation failed.")
                # Change back to original directory
                os.chdir(current_directory)
                raise SystemError(
                    os.system(
                        'tail -1 '
                        +os.path.join(
                            self.solver.case_folder,
                            'log.blockMesh')))
        elif self.solver.multiRegion:
            for region in self.solver.regions:
                code = os.system(
                    'blockMesh -region ' + region
                    + ' 2>&1 > log.blockMesh_'  +region)
                # Recreate coordinates
                code2 = os.system(
                    'postProcess -func writeCellCentres -region '
                    + region + ' 2>&1 > log.writeCellCentres_'+region)
                # Check if everything went well
                if (code == 0) & (code2 == 0):
                    print("Mesh for "+region+" region created!")
                elif (code != 0) & (code2 == 0):
                    print("ERROR: Mesh creation of "+region+" region failed.")
                    # Change back to original directory
                    os.chdir(current_directory)
                    raise SystemError(
                        os.system(
                            'tail -1 '
                            +os.path.join(
                                self.solver.case_folder,
                                'log.blockMesh_'+region)))
                elif (code == 0) & (code2 != 0):
                    print(
                        "ERROR: Coordinate creation of "
                        + region
                        + " region failed.")
                    # Change back to original directory
                    os.chdir(current_directory)
                    raise SystemError(
                        os.system(
                            'tail -1 '
                            +os.path.join(
                                self.solver.case_folder,
                                'log.writeCellCentres_'+region)))
                else:
                    print(
                        "ERROR: Both mesh and coordinate creation of "
                        + region
                        + " region failed.")
                    # Change back to original directory
                    os.chdir(current_directory)
                    raise SystemError(
                        os.system(
                            'tail -1 '
                            +os.path.join(
                                self.solver.case_folder,
                                'log.blockMesh_'+region)))
        # Change back to original directory
        os.chdir(current_directory)

    def writeCellCentres(self):
        """
        Writes the cell centres for all regions.
        """
        print("Writing cell centres ....")
        if not self.solver.multiRegion:
            code = os.system('postProcess -case ' + self.solver.case_folder +
                             ' -func writeCellCentres 2>&1 > log.writeCellCentres')
            if code!=0:
                print("ERROR: Writing cell centres failed")
                raise SystemError(os.system('tail -1 log.writeCellCentres'))
        elif self.solver.multiRegion:
            for region in self.solver.regions:
                code = os.system('postProcess -case '+self.solver.case_folder+' -func writeCellCentres -region ' +
                                 region + ' 2>&1 > log.writeCellCentres_' + region)
                if code != 0:
                    print("ERROR: Writing cell centres failed for region "+region)
                    raise SystemError(os.system('tail -1 log.writeCellCentres_'+region))
        print("Cell centres written")

    def getLatestTime(self):
        """
        Gets the latest time
        """
        dir_list = os.listdir(self.solver.case_folder)
        result_list = []
        for dir_name in dir_list:
            if dir_name.isdigit():
                result_list.append(float(dir_name))
        latest_folder_float = max(result_list)
        latest_folder = ('%f' % latest_folder_float).rstrip('0').rstrip('.')
        self.solver.latest_folder = latest_folder

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, indent=2)

    def setWriteFormat(self,write_format):
        """Writes the write format in the controlDict"""
        self.write_format = write_format
        if write_format == 'ascii' or 'binary':
            print("Setting write format to %s." %write_format)
            controlDict = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    'system',
                    'controlDict'))
            controlDict['writeFormat'] = write_format
            controlDict.writeFile()
            print("Write format set!")
        else:
            print("ERROR: Write format %s unknown." %write_format)
            raise

    def edit_decomposeParDict(self,region):
        """Edits the decomposeParDict and writes the cellDist"""
        # Edit decomposeParDict
        decomposeParDict = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'system',
                region,
                'decomposeParDict'
            )
        )
        decomposeParDict['numberOfSubdomains'] = self.solver.n_proc
        decomposeParDict['method'] = 'scotch'
        decomposeParDict.writeFile()
        # Write cellDist
        if region == '':
            # Single region case
            os.system(
                'decomposePar -cellDist > /dev/null 2>&1'
            )
        else:
            # Multi region case
            os.system(
                f'decomposePar -region {region} -cellDist '
                f'> /dev/null 2>&1'
            )
        os.system('rm -r processor* > /dev/null 2>&1')
        # Set decompose method back to manual
        decomposeParDict['method'] = 'manual'
        decomposeParDict.writeFile()

    def edit_setFieldsDict(self, region, direction='z'):
        """Edits the setFieldsDict"""
        # Edit setFieldsDict
        setFieldsDict = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'system',
                region,
                'setFieldsDict'
            ))
        small = 1e-8
        setFieldsDict['regions'] = []
        if direction == 'z':
            dz = self.geometry.h/self.solver.n_proc
        elif direction == 'x':
            blockMeshDict = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    'system',
                    region,
                    'blockMeshDict'
                ))
            x_coordinates = np.array(blockMeshDict['vertices'])[:,0]
            x_min = np.min(x_coordinates)
            x_max = np.max(x_coordinates)
            dx = (x_max - x_min)/self.solver.n_proc
        for i in range(self.solver.n_proc-1):
            setFieldsDict['regions'].append('boxToCell')
            reg = DictProxy()
            if direction == 'z':
                reg['box'] = TupleProxy(
                        ['(0 -1 %.8f)' %((dz)*(i+1)+small), 
                        '(1 1 %.8f)' %((dz)*(i+2)+small)]
                )
            elif direction == 'x':
                reg['box'] = TupleProxy(
                        ['(%.8f -1 0)' %(x_min+(dx)*(i+1)+small), 
                        '(%.8f 1 2)' %(x_min+(dx)*(i+2)+small)]
                )
            reg['fieldValues'] = ['volScalarFieldValue', 'cellDist', (i+1)]
            setFieldsDict['regions'].append(
                reg
            )
        setFieldsDict.writeFile()

        def repair_setFieldsDict():
            # Repair setFieldsDict by deleting the 21st line
            with open(os.path.join(
                    self.solver.case_folder,
                    'system',
                    region,
                    'setFieldsDict'
                    ), "r") as setFieldsDict:
                lines = setFieldsDict.readlines()
            del lines[20]
            with open(os.path.join(
                    self.solver.case_folder,
                    'system',
                    region,
                    'setFieldsDict'
                    ), "w") as setFieldsDict:
                setFieldsDict.writelines(lines)

        # Check if setFieldsDict is broken
        def check_setFieldsDict():
            try:
                setFieldsDict_check = ParsedParameterFile(
                    os.path.join(
                        self.solver.case_folder,
                        'system',
                        region,
                        'setFieldsDict'
                    ))
                return setFieldsDict_check
            except:
                print(
                    "ERROR: setFieldsDict broken. Please repair manually.")
                raise SystemError

        # Check setFieldsDict
        setFieldsDict_check = check_setFieldsDict()
        if type(setFieldsDict_check['regions'][0]) == float:
            print("ERROR: setFieldsDict not edited correctly.")
            print("Trying to repair ....")
            repair_setFieldsDict()
            # Check if repair worked
            check_setFieldsDict()
            print("Repair successful!")

    def write_cellAllocation(self,region):
        """Writes the cellDist file for manual decomposition"""
        # Read file
        cellDist_path = os.path.join(
            self.solver.case_folder,
            '0',
            region,
            'cellDist'
        )
        cellAllocation_path = os.path.join(
            self.solver.case_folder,
            'constant',
            region,
            'cellAllocation'
        )
        with open(cellDist_path, "r") as cellDist:
            lines = cellDist.readlines()
        # Check if all cells are in one processor
        uniform = False
        if lines[19].find("nonuniform") == -1:
            print(f"All cells in one processor.")
            uniform = True
            proc = int(lines[19].split()[-1][:-1])
            cellDecomposition_path = os.path.join(
                self.solver.case_folder,
                'constant',
                region,
                'cellDecomposition'
            )
            with open(cellDecomposition_path, "r") as cellDecomposition:
                lines_cellDecomposition = cellDecomposition.readlines()
            n_cells = int(lines_cellDecomposition[18].split()[0])

        # Replace and delete unwanted lines
        lines[11] = lines[11].replace("volScalarField", "labelList")
        lines[12] = lines[12].replace("0", "constant")
        lines[13] = lines[13].replace("cellDist", "cellAllocation")
        del lines[17:21]
        bfline = 0
        for i, line in enumerate(lines):
            if (line.find("boundaryField") != -1):
                bfline = i
                break
        del lines[bfline:]

        if uniform:
            # Append uniform distribution
            lines.append("(\n")
            for i in range(n_cells):
                lines.append(f"{proc}\n")
            lines.append(")\n")
            lines.append(";\n")
            lines.append("\n")

        # Write file
        with open(cellAllocation_path, "w") as cellAllocation:
            for line in lines:
                cellAllocation.write(line)

    def decomposeCase(self, direction='z'):
        """Decomposes the case"""
        if self.solver.n_proc == 1:
            print("ERROR: n_proc = 1. No decomposition necessary.")
            return
        assert direction in ['x', 'z'], "ERROR: direction must be 'x' or 'z'"
        print("Decomposing case ....")
        # Set write format to ascii
        original_write_format = self.write_format
        self.setWriteFormat('ascii')
        # Save current working directory to go back afterwards
        current_directory = os.getcwd()
        # Got to case folder
        os.chdir(self.solver.case_folder)
        os.system('rm -r processor* > /dev/null 2>&1')

        for region in self.solver.regions:
            print(f"Preparing decomposition of region {region} ....")
            # Prepare decomposition
            self.edit_setFieldsDict(region, direction=direction)            
            self.edit_decomposeParDict(region)
            # Create setFields file
            os.system(
                f'setFields -region {region} > log.setFields_{region} 2>&1'
            )
            print(f"Done!")
            print("Writing cell allocation for  region %s ...." %region)
            # Write cellAllocation file
            self.write_cellAllocation(region)
            print("Cell allocation of region %s written!" %region)
        if not self.solver.multiRegion:
            print(f"Preparing decomposition ....")
            # Prepare decomposition
            self.edit_setFieldsDict('', direction=direction)
            self.edit_decomposeParDict('')
            # Create setFields file
            os.system(
                'setFields > log.setFields 2>&1'
            )
            print(f"Done!")
            print("Writing cell allocation....")
            # Write cellAllocation file
            self.write_cellAllocation('')
            print("Cell allocation written!")

        self.setWriteFormat(original_write_format)
        for region in self.solver.regions:
            print("Decomposing region %s ...." %region)
            # Decompose case
            os.system(
                f'decomposePar -region {region} '
                f'> log.decomposePar_{region} 2>&1'
            )
            print("Decomposition of region %s done!" %region)
        if not self.solver.multiRegion:
            print("Decomposing....")
            # Decompose case
            os.system(
                'decomposePar > log.decomposePar 2>&1'
            )
            print("Decomposition done!")

        for proc in range(self.solver.n_proc):
            os.system(
                f'touch processor{proc}/'
                f'{self.solver.working_folder_name}_{proc}.foam'
            )

        # Change back to original directory
        os.chdir(current_directory)

    def setMaterialProperties(self):
        """
        Sets the material properties as defined in self.material
        """
        print("Setting material properties ....")
        # ----------- Air side ----------- #
        print("Setting air properties ....")
        airProp = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'constant',
                'air',
                'thermophysicalProperties'))
        airProp['mixture']['equationOfState']['rho'][2] = (
            self.material.air.rho)
        airProp['mixture']['transport']['lambda'][2] = (
            self.material.air.kappa)
        airProp['mixture']['transport']['TDiff'][2] = (
            self.material.air.alpha)
        airProp['mixture']['transport']['cwDiff'][2] = (
            self.material.air.cwDiff)
        airProp['mixture']['transport']['p0'][2] = (
            self.material.air.p0)
        airProp['mixture']['transport']['beta'][2] = (
            self.material.air.beta)
        airProp['mixture']['transport']['nu'][2] = (
            self.material.air.nu)
        airProp['mixture']['humidity']['deltaH'][2] = (
            self.material.air.deltaH)
        airProp.writeFile()
        try:
            controlDict = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    'system',
                    'controlDict'))
            controlDict['functions']['totalPressure']['rhoInf'] = (
                self.material.air.rho)
            controlDict.writeFile()
        except:
            print("No air density in controlDict.")
        print("Air properties set!")
        # ----------- Solution side ----------- #
        for region in self.solver.solutionRegions:
            print("Setting %s properties ...." %region)
            solProp = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    'constant',
                    region,
                    'thermophysicalProperties'))
            # Set material
            if self.material.solution.medium == 'LiBr':
                # Use Patterson (1984) correlation
                solProp['mixture']['transport']['flagMedium'] = 0
            elif self.material.solution.medium == 'custom':
                # Use custom Antoine fit
                solProp['mixture']['transport']['flagMedium'] = 1
            else:
                print("ERROR: No known solution media set.")  
                raise
            # Set thermophysical data
            solProp['mixture']['equationOfState']['rho'][2] = (
                self.material.solution.rho)
            solProp['mixture']['transport']['lambda'][2] = (
                self.material.solution.kappa)
            solProp['mixture']['transport']['TDiff'][2] = (
                self.material.solution.alpha)
            solProp['mixture']['transport']['cwDiff'][2] = (
                self.material.solution.cwDiff)
            solProp['mixture']['transport']['p0'][2] = (
                self.material.solution.p0)
            solProp['mixture']['transport']['nu'][2] = (
                self.material.solution.nu)
            # Set fit (if flagFit=True)
            if self.material.solution.flagFit:
                self.setSolutionFitParameters(
                    solutionPropertiesDict=solProp,
                    region=region)
            solProp.writeFile()
            print("%s properties set!" %region)
        print("Material properties set!")

    def setSolutionFitParameters(self, 
            solutionPropertiesDict=False, region='solution'):
        """
        Sets the solution fit properties as defined in 
        self.material.solution.fit
        """
        print("Setting solution fit parameters ....")
        if not solutionPropertiesDict:
            solProp = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    'constant',
                    region,
                    'thermophysicalProperties'))
        else:
            solProp = solutionPropertiesDict
        solProp['mixture']['vaporPressureFit']['k0A'][2] = (
            self.material.solution.fit.k0A)
        solProp['mixture']['vaporPressureFit']['k1A'][2] = (
            self.material.solution.fit.k1A)
        solProp['mixture']['vaporPressureFit']['k2A'][2] = (
            self.material.solution.fit.k2A)
        solProp['mixture']['vaporPressureFit']['k3A'][2] = (
            self.material.solution.fit.k3A)
        solProp['mixture']['vaporPressureFit']['k0B'][2] = (
            self.material.solution.fit.k0B)
        solProp['mixture']['vaporPressureFit']['k1B'][2] = (
            self.material.solution.fit.k1B)
        solProp['mixture']['vaporPressureFit']['k2B'][2] = (
            self.material.solution.fit.k2B)
        solProp['mixture']['vaporPressureFit']['k3B'][2] = (
            self.material.solution.fit.k3B)
        solProp['mixture']['vaporPressureFit']['k0C'][2] = (
            self.material.solution.fit.k0C)
        solProp['mixture']['vaporPressureFit']['k1C'][2] = (
            self.material.solution.fit.k1C)
        solProp['mixture']['vaporPressureFit']['k2C'][2] = (
            self.material.solution.fit.k2C)
        solProp['mixture']['vaporPressureFit']['k3C'][2] = (
            self.material.solution.fit.k3C)
        solProp['mixture']['vaporPressureFit']['D'][2] = (
            self.material.solution.fit.D)
        solProp.writeFile()
        print("Solution parameters set!")

    def setAdiabatic(self):
        """
        Sets the temperature boundary condition adiabatic. Call function after
        self.setInitialization()!
        """
        print("Setting adiabatic boundary condition ....")
        for region in self.solver.solutionRegions:
            T_sol = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    str(self.solver.startTime),
                    region,
                    'T'))
            T_sol['boundaryField']['wall']['type'] = 'zeroGradient'
            T_sol['boundaryField']['wall']['value'] = ('uniform '
                + str(self.initialization.solution.T_wall+273.15))
            T_sol.writeFile()
        print("Adiabatic boundary condition set!")

    def setDesorption(self):
        """
        Sets the vapor pressure to desorption. Call function after
        self.setInitialization()!
        """
        print("Setting desorption conditions ....")
        for region in self.solver.solutionRegions:
            pw_sol = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    str(self.solver.startTime),
                    region,
                    'pw'))
            pw_sol['internalField'] = ('uniform '
                + str(self.initialization.air.pw_in*2))
            pw_sol.writeFile()
        print("Desorption conditions set!")

    def setNoMassTransfer(self):
        """
        Sets the solver to no mass transfer. Call function after
        self.setInitialization()!
        """
        print("Setting no mass transfer conditions ....")
        # Set cw in air region
        C_air = ParsedParameterFile(
            os.path.join(self.solver.case_folder,
            str(self.solver.startTime),
            'air', 
            'cw'))
        for sol_reg in self.solver.solutionRegions:
            C_air['boundaryField'][f'air_to_{sol_reg}']['type'] = (
                'fixedValue')
            C_air['boundaryField'][f'air_to_{sol_reg}']['value'] = (
                '$internalField')
            del (C_air['boundaryField']
                    [f'air_to_{sol_reg}']
                    ['neighbourRegionName'])
            del (C_air['boundaryField']
                    [f'air_to_{sol_reg}']
                    ['neighbourPatchName'])
            del (C_air['boundaryField']
                    [f'air_to_{sol_reg}']
                    ['neighbourFieldName'])
            del (C_air['boundaryField']
                    [f'air_to_{sol_reg}']
                    ['relaxationFactor'])
        C_air.writeFile()
        # Set cw in solution region
        for region in self.solver.solutionRegions:
            C_sol = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    str(self.solver.startTime),
                    region,
                    'cw'))
            C_sol['boundaryField']['solution_to_air']['type'] = (
                'fixedValue')
            C_sol['boundaryField']['solution_to_air']['value'] = (
                '$internalField')
            del C_sol['boundaryField']['solution_to_air']['neighbourRegionName']
            del C_sol['boundaryField']['solution_to_air']['neighbourPatchName']
            del C_sol['boundaryField']['solution_to_air']['neighbourFieldName']
            C_sol.writeFile()
        print("No mass transfer conditions set!")

    def makeSingleRegion(self):
        """
        Makes the case single region
        """
        print("Making the case single region ....")
        # 0 folder
        # Moving the old folder
        os.system(
            'mv '
            + os.path.join(
                self.solver.case_folder,
                '0'
            )
            + ' '
            + os.path.join(
                self.solver.case_folder,
                '0_multi'
            )
        )
        # Creating the new 0 folder
        os.system(
            'mkdir '
            + os.path.join(
                self.solver.case_folder,
                '0'
            )
        )
        # Copy U and p
        os.system(
            'cp '
            + os.path.join(
                self.solver.case_folder,
                '0_multi',
                'air',
                'U'
            )
            + ' '
            + os.path.join(
                self.solver.case_folder,
                '0'
            )
        )
        os.system(
            'cp '
            + os.path.join(
                self.solver.case_folder,
                '0_multi',
                'air',
                'p'
            )
            + ' '
            + os.path.join(
                self.solver.case_folder,
                '0'
            )
        )
        # Remove old 0 folder
        os.system(
            'rm -r '
            + os.path.join(
                self.solver.case_folder,
                '0_multi'
            )
        )
        # constant folder
        # Moving the old folder
        os.system(
            'mv '
            + os.path.join(
                self.solver.case_folder,
                'constant'
            )
            + ' '
            + os.path.join(
                self.solver.case_folder,
                'constant_multi'
            )
        )
        # Creating the new constant folder
        os.system(
            'mkdir '
            + os.path.join(
                self.solver.case_folder,
                'constant'
            )
        )
        # Copying the old constant folder
        os.system(
            'cp -r '
            + os.path.join(
                self.solver.case_folder,
                'constant_multi',
                'air',
                '*'
            )
            + ' '
            + os.path.join(
                self.solver.case_folder,
                'constant'
            )
        )
        # Changing nu
        transportProperties = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'constant',
                'transportProperties'))
        transportProperties['nu'][1] = self.material.air.nu
        transportProperties.writeFile()
        # Removing the old constant folder
        os.system(
            'rm -r '
            + os.path.join(
                self.solver.case_folder,
                'constant_multi'
            )
        )
        # system folder
        # Moving the old folder
        os.system(
            'mv '
            + os.path.join(
                self.solver.case_folder,
                'system'
            )
            + ' '
            + os.path.join(
                self.solver.case_folder,
                'system_multi'
            )
        )
        # Creating the new system folder
        os.system(
            'mkdir '
            + os.path.join(
                self.solver.case_folder,
                'system'
            )
        )
        # Copying the old system folder
        os.system(
            'cp -r '
            + os.path.join(
                self.solver.case_folder,
                'system_multi',
                'air',
                '*'
            )
            + ' '
            + os.path.join(
                self.solver.case_folder,
                'system'
            )
        )
        os.system(
            'cp -r '
            + os.path.join(
                self.solver.case_folder,
                'system_multi',
                'controlDict'
            )
            + ' '
            + os.path.join(
                self.solver.case_folder,
                'system'
            )
        )
        # Removing the old system folder
        os.system(
            'rm -r '
            + os.path.join(
                self.solver.case_folder,
                'system_multi'
            )
        )
        # Editing the controlDict
        controlDict = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'system',
                'controlDict'))
        if self.classtype == 'wire':
            controlDict['application'] = 'pimpleFoam'
        elif self.classtype == 'stackedArrangement':
            controlDict['application'] = 'simpleFoam'
        del controlDict['functions']['residuals_air']
        for region in self.solver.solutionRegions:
            del controlDict['functions'][f'residuals_{region}']
        del controlDict['functions']['residuals_flow_air']['region']
        del controlDict['functions']['absorbedMassFluxAir']
        del controlDict['functions']['absorbedMassFluxAirObject']
        del controlDict['functions']['sensibleHeatFluxAir']
        del controlDict['functions']['sensibleHeatFluxDensityAirObject']
        del controlDict['functions']['outletTemperature']
        del controlDict['functions']['staticPressureDifference']['region']
        del controlDict['functions']['totalPressure']['region']
        del controlDict['functions']['totalPressureDifference']['region']
        controlDict.writeFile()
        # Edit fvSolutions
        fvSolutions = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'system',
                'fvSolution'))
        fvSolutions['relaxationFactors']['fields']['p'] = 0.3
        fvSolutions['relaxationFactors']['equations']['U'] = 0.7
        fvSolutions.writeFile()
        # Edit fvSchemes
        if self.classtype == 'stackedArrangement':
            fvSchemes = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    'system',
                    'fvSchemes'))
            fvSchemes['ddtSchemes']['default'] = 'steadyState'
            fvSchemes.writeFile()
        # Redo decomposition
        # Remove old processor folders
        os.system(
            'rm -r '
            + os.path.join(
                self.solver.case_folder,
                'processor*'
            )
        )
        # Redo decomposition
        os.system(
            'decomposePar -case '
            + self.solver.case_folder
            + ' > '
            + os.path.join(
                self.solver.case_folder,
                'log.decomposePar'
            )
            + ' 2>&1'
        )
        print("Case made single region!")

    def makeSteadyState(self):
        """
        Make the case steady state.
        """
        # Change controlDict
        controlDict = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'system',
                'controlDict'))
        controlDict['application'] = 'heatMassTransferFoam'
        # Delete functions
        del controlDict['functions']['residuals_flow_air']
        del controlDict['functions']['staticPressureDifference']
        del controlDict['functions']['totalPressure']
        del controlDict['functions']['totalPressureDifference']
        controlDict.writeFile()
        # Change fvSolution
        for region in self.solver.regions:
            fvSolution = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    'system',
                    region,
                    'fvSolution'))
            del fvSolution['HMT']['residualControl']['T']
            del fvSolution['HMT']['residualControl']['cw']
            fvSolution.writeFile()
        # Change fvSchemes
        for region in self.solver.regions:
            fvSchemes = ParsedParameterFile(
                os.path.join(
                    self.solver.case_folder,
                    'system',
                    region,
                    'fvSchemes'))
            fvSchemes['ddtSchemes']['default'] = 'steadyState'
            fvSchemes.writeFile()

    def setMaxCo(self, maxCo = 0.7):
        """
        Set maximum Courant number in controlDict.
        """
        # Change controlDict
        controlDict = ParsedParameterFile(
            os.path.join(
                self.solver.case_folder,
                'system',
                'controlDict'))
        # Set max Co
        controlDict['adjustTimeStep'] = 'yes'
        controlDict['maxCo'] = maxCo
        controlDict['writeControl'] = 'adjustableRunTime'
        controlDict.writeFile()


if __name__ == "__main__":
    pass
