#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 08:36:28 2023

@author: cghiaus

Obtain the .csv file for inputs from 03CubicBuilding.py
This inputs file is to be used to know what to obtain from weather file and
information on walls the solar radiation on walls


For step input:
    in 03CubicBuilding.py:
u = np.zeros([8, n])                # u = [To To To Tisp Φo Φi Qa Φa]
u[0:3, :] = 10 * np.ones([3, n])    # To = 10 for n time steps
u[3, :] = 20 * np.ones([1, n])      # Tisp = 20 for n time steps

pd.DataFrame(u).to_csv('u_step.csv', index=False)
...
    in this file
u_read = pd.read_csv('u_step.csv')

For weather data:
    in 03CubicBuilding.py:
data.rename(columns={'temp_air': 'To'}).to_csv('u_weather.csv', index=False)
...

For input u for weather data:
    in 03CubicBuilding.py:

# input vector
To = data['To']
Ti = data['Ti']
Φo = α_wSW * wall['Surface']['Layer_out'] * data['Φtot']
Φi = τ_gSW * α_wSW * wall['Surface']['Glass'] * data['Φtot']
Qa = data['Qa']
Φa = α_gSW * wall['Surface']['Glass'] * data['Φtot']

u = pd.concat([To, To, To, Ti, Φo, Φi, Qa, Φa], axis=1)
u.columns.values[[4, 5, 7]] = ['Φo', 'Φi', 'Φa']

pd.DataFrame(u).to_csv('u_weather.csv', index=False)
pd.DataFrame(u).to_csv('u_weather.csv', index=False)
...
"""
import pandas as pd

file_path = 't03CubicBuilding.py'
# Read the contents of the script_to_run.py file
with open(file_path, 'r') as file:
    script_contents = file.read()

# Execute the script
exec(script_contents)


# Run the Python file using subprocess
u = pd.read_csv('u_step.csv')
