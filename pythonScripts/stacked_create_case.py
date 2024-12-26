#!/usr/bin/env python3
"""
01.06.2023
Felix Hochwallner
Calculate stacked arrangement simulations
"""

# Import standard python modules
import os
import numpy as np
import pandas as pd

# Import custom modules
from caseObjects.stackedArrangement import stackedArrangement

from caseObjects.settings.solver import solverSettings
from caseObjects.settings.geometry import (
    airStackedArrangementGeometrySettings,
    solutionStackedArrangementGeometrySettings,
    stackedArrangementGeometrySettings)
from caseObjects.settings.material import (
    airMaterialSettings,
    solutionMaterialSettings, 
    airSolutionMaterialSettings)
from caseObjects.settings.initialization import (
    airInitializationSettings,
    solutionInitializationSettings,
    airSolutionInitializationSettings)

from caseObjects.fluidProperties import humidAirProperties
from caseObjects.fluidProperties import waterProperties
from caseObjects.fluidProperties import LiBrProperties as desiccantProperties

# ----- Inputs -----
# Air inlet conditions
u_a_in = 1.3 # m/s
p_a_in = 1e5 # Pa
T_a_in = 4 # degC
RH_a_in = 80/100 # 1

# Solution inlet conditions
desiccant = 'LiBr'
T_s_in = 2 # degC
T_wall = 2 # degC
Chi_s_in = 0.5 # 1
delta_Nu = 0.75e-3 # m

# Geometry
l = 80e-3 # m /Air gap length
b = 6e-3 # m / Air gap width
h = 500e-3 # m / Air gap height
n_steps = 2 # 1 / Number of steps
d_wall = 0.4e-3 # m / Wall thickness
grading = 10 # 1 / Grading of mesh in solution film
x_gap = 3e-3 # m / Gap between steps
x_inlet = 1.5 * b # m / Inlet length
x_outlet = 30 * b # m / Outlet length

# Mesh
overwrite_check = True
calc_ny_sol = False
ny_air = 35 # 1 / Perpendicular to length
ny_sol = 10 # 1 / Along film thickness
nz = 1 # 1 / Height direction
grading_inlet_outlet = True
grading_inlet_outlet_ratio = 1.1
n_proc = 4 # 1 / Number of processes

# Solver
t_end = 10000 # s
time_step = 1 # s
write_interval = 10000 # s

# ----- Solver settings -----
case_folder = os.path.abspath(
    os.path.join(
        '..',
        'baseCases',
        'stacked'))
solver = 'simpleFoam'
results_folder = ''
working_folder = os.path.abspath(
    os.path.join(
        'workingFolder'))
# Check if working folder exists and create it if not
if not os.path.exists(working_folder):
    os.makedirs(working_folder)

solver_settings = solverSettings(
    case_folder, 
    solver, 
    results_folder,
    working_folder = working_folder)
solver_settings.setTimes(
    endTime = t_end,
    timeStep = time_step,
    writeInterval = write_interval)

# ----- Geometry settings -----
l_step = (l - x_gap * (n_steps-1))/n_steps
air_geometry = airStackedArrangementGeometrySettings(
    l=l_step, 
    b=b, 
    h=h,
    d_wall=d_wall,
    n_steps=n_steps,
    ny_air=ny_air,
    ny_sol=ny_sol,
    nz=nz,
    x_gap=x_gap,
    x_inlet=x_inlet,
    x_outlet=x_outlet,
    grading_inlet_outlet=grading_inlet_outlet,
    grading_inlet_outlet_ratio=grading_inlet_outlet_ratio)

solution_geometry = solutionStackedArrangementGeometrySettings(
    l=l_step, 
    b=b, 
    h=h,
    d_wall=d_wall,
    n_steps=n_steps,
    ny_air=ny_air,
    ny_sol=ny_sol,
    nz=nz,
    grading=grading,
    x_gap=x_gap,
    x_inlet=x_inlet,
    x_outlet=x_outlet,
    grading_inlet_outlet=grading_inlet_outlet,
    grading_inlet_outlet_ratio=grading_inlet_outlet_ratio)

geometry_settings = stackedArrangementGeometrySettings(
    air_geometry, 
    solution_geometry,
    overwrite_check=overwrite_check,
    calc_ny_sol=calc_ny_sol)

# ----- Material settings ------
# Air
air_material = airMaterialSettings()
air_material.rho = humidAirProperties.rho(p_a_in, T_a_in, RH_a_in)
air_material.cp = humidAirProperties.cp(p_a_in, T_a_in, RH_a_in)
air_material.mu = humidAirProperties.mu(p_a_in, T_a_in, RH_a_in)
air_material.Pr = humidAirProperties.Pr(p_a_in, T_a_in, RH_a_in)
air_material.kappa = humidAirProperties.kappa(p_a_in, T_a_in, RH_a_in)
air_material.cwDiff = humidAirProperties.cwDiff(
    p_a_in, T_a_in, which='VDI')
air_material.p0 = p_a_in # Pa
air_material.beta = humidAirProperties.beta()
air_material.deltaH = waterProperties.hVap((T_s_in+T_a_in)/2)
air_material.nu = humidAirProperties.nu(p_a_in, T_a_in, RH_a_in)
air_material.alpha = humidAirProperties.alpha(p_a_in, T_a_in, RH_a_in)
# Solution
solution_material = solutionMaterialSettings(desiccant)
solution_material.medium = desiccant
solution_material.rho = desiccantProperties.rho(Chi_s_in, T_s_in)
solution_material.cp = desiccantProperties.cp(Chi_s_in, T_s_in)
solution_material.mu = desiccantProperties.mu(Chi_s_in, T_s_in)
solution_material.kappa = desiccantProperties.kappa(Chi_s_in, T_s_in)
solution_material.cwDiff = desiccantProperties.cwDiff(Chi_s_in, T_s_in)
solution_material.p0 = p_a_in # Pa
solution_material.hs = waterProperties.hVap((T_s_in+T_a_in)/2)
solution_material.nu = desiccantProperties.nu(Chi_s_in, T_s_in)
solution_material.alpha = desiccantProperties.alpha(Chi_s_in, T_s_in)
solution_material.setFit()
material_settings = airSolutionMaterialSettings(
    air_material, 
    solution_material)

# ----- Initialization -----
C_a_in = humidAirProperties.cwRH(
    p_a_in, 
    T_a_in, 
    RH_a_in)
C_s_in = 1 - Chi_s_in

air_initialization = airInitializationSettings(
    v_in = u_a_in,
    T_in = T_a_in,
    C_in = C_a_in
)
rho_s = solution_material.rho
nu_s = solution_material.nu
g = 9.81 # m/s^2
m_dot_s = (l*rho_s*g)/(3*nu_s)*delta_Nu**3
solution_initialization = solutionInitializationSettings(
    T_in = T_s_in, 
    C_in = C_s_in, 
    T_wall = T_wall, 
    mDot = m_dot_s
)
initialization_settings = airSolutionInitializationSettings(
    air_initialization,
    solution_initialization
)

# ----- Define case -----
case = 'stacked'
stacked_case = stackedArrangement(
    solver_settings,
    geometry_settings,
    material_settings,
    initialization_settings,
    case=case
)
stacked_case.setBlockMeshDict()
stacked_case.setWriteFormat(write_format='binary')
stacked_case.createMesh()
stacked_case.setMaterialProperties()
stacked_case.setInitialization()
stacked_case.solver.n_proc = n_proc
if n_proc > 1:
    stacked_case.decomposeCase(direction='x')
stacked_case.makeSteadyState()
