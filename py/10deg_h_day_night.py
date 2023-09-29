#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 14:47:05 2023

@author: cghiaus
"""
import pandas as pd
from dm4bem import read_epw

# Inputs
# ======
filename = './weather_data/FRA_Lyon.074810_IWEC.epw'

θ = 18          # °C, indoor temperature all time
θday = θ        # °C, indoor temperature during day: 06:00 - 23:00
θnight = 16     # °C, indoor temperature during night 23:00 - 06:00

period_start = '2000-01-01'
period_end = '2000-12-31'

daytime_start = '06:00:00+01:00'
daytime_end = '22:00:00+01:00'

# Computation
# ===========
# read Energy Plus Weather data (file .EPW)
[data, meta] = read_epw(filename, coerce_year=2000)

# select outdoor air temperature; call it θout
df = data[["temp_air"]]
del data
df = df.rename(columns={'temp_air': 'θout'})

# Select the data for a period of the year
df = df.loc[period_start:period_end]

# Compute degree-hours for fixed set-point
# ----------------------------------------
df['Δθfix'] = θ - df['θout'].where(θ > df['θout'], θ)

# Compute degree-hours for variable (day/night) set-point
# -------------------------------------------------------
# Define start time for day and night
day_start = pd.to_datetime(daytime_start).time()
day_end = pd.to_datetime(daytime_end).time()

# Daytime should be between 00:00 and 24:00
# Daytime including midnight is not allowed, e.g., 22:00 till 06:00
day = (df.index.time >= day_start) & (df.index.time <= day_end)
night = ~day

# Degree-hours for daytime
df['Δθday'] = θday - df['θout'].where(
    (θday > df['θout']) & day,
    θday)

# Degree-hours for nighttime
df['Δθnight'] = θnight - df['θout'].where(
    (θnight > df['θout']) & night,
    θnight)

# Sum of degree-hours for fix indoor temperature
dh_fix = df['Δθfix'].sum()

# Sum of degree-hours for intermitent heating
dh_interm = df['Δθday'].sum() + df['Δθnight'].sum()

# Results
# =======
print(f"degree-hours fix set-point: {dh_fix:.1f}")
print(f"degree-hours variable set-point: {dh_interm:.1f}")
print(f"Estimated savings: {(dh_fix - dh_interm) / dh_interm * 100:.0f} %")
