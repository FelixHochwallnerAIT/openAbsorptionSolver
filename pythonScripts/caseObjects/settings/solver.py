#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:17:21 2020

@author: Felix Hochwallner

Class definition solverSettings, for OpenFoam case classes.
"""
import os
from datetime import datetime
import numpy as np
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

class solverSettings:
    """
    Stores the solver settings of an OpenFoam case.
    """
    def __init__(self, 
            folder, 
            solver, 
            results, 
            working_folder = '.',
            transient = False,
            n_proc = 1,
            n_nodes = 1,
            startTime = 0.0,
            endTime = 100.0,
            timeStep = 1.0,
            writeInterval = 1,
            cluster = 'None'):
        """
        Constructor class
        Input:
        folder: Case folder location
        solver: OpenFoam solver used
        results: Results folder name
        """
        self.case_folder = os.path.abspath(folder)
        now = datetime.now()
        current_time = now.strftime("%Y-%m-%d_%H:%M:%S")
        workingFolderName = (
            os.path.basename(
                os.path.normpath(self.case_folder))+
                '_'+current_time)
        self.working_folder = os.path.abspath(
            os.path.join(
                '.',
                workingFolderName))
        # Make path absolute and create folder, copy content ...
        self.solver = solver
        self.results_folder = results
        self.working_folder = working_folder
        self.latest_folder = str(startTime)
        self.transient = transient
        self.n_proc = n_proc
        self.n_nodes = n_nodes
        self.startTime = startTime
        self.endTime = endTime
        self.timeStep = timeStep
        self.writeInterval = writeInterval
        self.cluster = cluster
        # Check if case is multi region
        try:
            regionProperties = ParsedParameterFile(
                os.path.join(
                    self.case_folder, 
                    'constant', 
                    'regionProperties')
            )
            regions = regionProperties['regions'][1]
            print("Multiple regions found: %s" %regions)
            self.multiRegion = True
            self.regions = regions
            self.solutionRegions = (
                [reg for reg in regions if 'solution' in reg])
            self.airRegions = (
                [reg for reg in regions if 'air' in reg])
        except:
            print("Single region case found.")
            self.multiRegion = False
            self.regions = []

    def setTimes(self, startTime=0, endTime=0, timeStep=1, writeInterval=1):
        """
        Set solver times.
        Input:
        startTime
        endTime
        timeStep
        writeInterval
        """
        self.startTime = startTime
        self.endTime = endTime
        self.timeStep = timeStep
        self.writeInterval = writeInterval
        self.setControlDict()

    def setTransient(self):
        """
        Sets solver transient
        """
        self.transient = True

    def setSteadyState(self):
        """
        Sets solver steady state
        """
        self.transient = False


    def setCluster(self, cluster):
        """
        Sets the cluster
        """
        self.cluster = cluster

    def setControlDict(self):
        """
        Setting controlDict.
        """
        print("Setting controlDict.")
        controlDict = ParsedParameterFile(self.case_folder + '/system/controlDict')
        controlDict['startTime'] = self.startTime
        controlDict['deltaT'] = self.timeStep
        controlDict['endTime'] = self.endTime
        controlDict['writeInterval'] = self.writeInterval
        if self.writeInterval < 1:
            controlDict['writeControl'] = 'adjustableRunTime'
        controlDict.writeFile()
        print("controlDict set!")
