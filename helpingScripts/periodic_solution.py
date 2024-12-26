#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
27.02.2024
Felix Hochwallner
Extract periodic solution
"""

#%%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from Ofpp import parse_internal_field, parse_field_all, parse_boundary_field
import json

# %%

path_data = './workingFolder/'
folder = 'wire_init'

# Get time steps
time_steps = os.listdir(
    os.path.join(
        path_data,
        folder)
    )
# Remove non-numeric entries
time_steps = [x for x in time_steps if x.replace('.', '', 1).isdigit()]

# Remove 0 time step
time_steps.remove('0')

# Remove all time steps before 0.1
time_steps = [x for x in time_steps if float(x) >= 0.1]

# Sort time steps
time_steps = sorted(time_steps, key=float)

# %%
# Get pressure fields
p_fields = dict()
for time_step in time_steps:
    print("\nReading time step: ", time_step, " ....")
    p_field = parse_internal_field(
        os.path.join(
            path_data,
            folder,
            time_step,
            'p'))
    p_fields[float(time_step)] = p_field

# %%
# Get pressure difference fields
p_diff = dict()
p_diff_mean = dict()
for time_step in time_steps:
    print("\nCalculating pressure difference at time step: ", 
        time_step, 
        " ....")
    p_diff[float(time_step)] = (
        p_fields[float(time_step)]
        - p_fields[float(time_steps[0])])
    p_diff_mean[float(time_step)] = np.mean(p_diff[float(time_step)])

# %%
# Find 0 crossings of p_diff_mean
p_diff_mean_crossings = np.where(
    np.diff(np.sign(list(p_diff_mean.values()))))[0]


# %%
# Plot pressure difference
plt.figure()
plt.plot(
    p_diff_mean.keys(),
    p_diff_mean.values())
plt.scatter(
    [list(p_diff_mean.keys())[x] for x in p_diff_mean_crossings],
    [list(p_diff_mean.values())[x] for x in p_diff_mean_crossings],
    color='red')
plt.xlabel('Time (s)')
plt.ylabel('Difference to 0 field (Pa)')
plt.grid()

# %%
# Take data until fourth crossing and repeat
time_steps_repeat = time_steps[:p_diff_mean_crossings[4]]
# Repeat a few times
# Create new time steps
delta_t = float(time_steps_repeat[-1])-float(time_steps_repeat[0])
time_steps_corresponding = time_steps[:p_diff_mean_crossings[4]]
for i in range(5):
    time_steps_corresponding += time_steps[:p_diff_mean_crossings[4]]
    time_steps_repeat += [
        str(round(float(x) + delta_t*(i+1), 8)) 
            for x in time_steps[:p_diff_mean_crossings[4]]]
    
# Create repeated pressure fields
p_fields_repeat = dict()
p_diff_mean_repeat = dict()
for i in range(len(time_steps_repeat)):
    p_fields_repeat[float(time_steps_repeat[i])] = (
        p_fields[float(time_steps_corresponding[i])])
    p_diff_mean_repeat[float(time_steps_repeat[i])] = (
        p_diff_mean[float(time_steps_corresponding[i])])

# %%
# Plot pressure difference vs repeated pressure difference
plt.figure()
plt.plot(
    p_diff_mean_repeat.keys(),
    p_diff_mean_repeat.values(),
    label='Repeated')
plt.plot(
    p_diff_mean.keys(),
    p_diff_mean.values(),
    label='Original')
plt.xlabel('Time (s)')
plt.ylabel('Difference to 0 field (Pa)')
plt.legend()
plt.grid()

# %%
# Remove all time steps above 4th crossing
time_steps_remove = time_steps[p_diff_mean_crossings[4]:]
for time_step in time_steps_remove:
    print("\nRemoving time step: ", time_step, " ....")
    os.system(
        'rm -r ' 
        + os.path.join(
            path_data,
            folder,
            time_step))

# %%

# 