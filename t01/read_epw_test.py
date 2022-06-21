#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 10:48:03 2021

@author: cghiaus
"""

import numpy as np
from dm4bem import read_epw

filename = 'FRA_Lyon.074810_IWEC.epw'
[data, meta] = read_epw(filename, coerce_year=None)
print(f"Years in the dataset: {sorted(set(data['year']))}")

weather_data = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
weather_data = weather_data[(weather_data.index >= '1983-1-1') & (
    weather_data.index < '1983-1-6')]
data

surface_orientation = {'slope': 90,
                       'azimuth': 0,
                       'latitude': 45}
albedo =  0.2

B = surface_orientation['slope']
Z = surface_orientation['azimuth']
L = surface_orientation['latitude']

# Transform degrees in radiants
B = B * np.pi / 180
Z = Z * np.pi / 180
L = L * np.pi / 180

n = weather_data.index.dayofyear

# eq. 1.6.1a
declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
d = declination_angle * np.pi / 180

# Example 1.6.1
hour = weather_data.index.hour
minute = weather_data.index.minute + 60
hour_angle = ((hour + minute / 60) - 12) * 15
h = hour_angle * np.pi / 180

# incidence angle eq. 1.6.2
theta = np.sin(d) * np.sin(L) * np.cos(B)
theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
theta -= np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
theta = np.array(np.arccos(theta))
theta[theta > (np.pi / 2)] = np.pi / 2

# Th-CE 2005
dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
dir_rad[dir_rad < 0] = 0

# Th-CE 2005 Eq. 79
dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(theta)) / 2

gamma = np.cos(d) * np.cos(L) * np.cos(h)
gamma += np.sin(d) * np.sin(L)
gamma = np.array(np.arcsin(gamma))
gamma[gamma < 1e-4] = 1e-4

# direct radiation on horizontal
dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)

# radiation reflected by the ground
ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
ref_rad *= (1 - np.cos(B) / 2)
