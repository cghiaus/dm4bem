#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 08:32:19 2021

@author: cghiaus
"""

from dm4bem import read_epw, sol_rad_tilt_surf

start_date = '2000-04-10'
end_date = '2000-05-16'

# Read weather data from Energyplus file
filename = 'FRA_Lyon.074810_IWEC.epw'
filename = 'FRA_AR_Lyon-Bron.AP.074800_TMYx.2004-2018.epw'

[data, meta] = read_epw(filename, coerce_year=None)
month_year = data['month'].astype(str) + '-' + data['year'].astype(str)
print(f"Months - years in the dataset: {sorted(set(month_year))}")

weather_data = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
weather_data.index = weather_data.index.map(lambda t: t.replace(year=2000))

weather_data = weather_data[(weather_data.index >= start_date) & (
    weather_data.index < end_date)]

# del data

weather_data.plot()

# Solar radiation on a tilted surface
surface_orientation = {'slope': 0,
                       'azimuth': 0,
                       'latitude': 45}
albedo = 0.2

rad_surf = sol_rad_tilt_surf(
    weather_data, surface_orientation, albedo)

rad_surf.plot()
