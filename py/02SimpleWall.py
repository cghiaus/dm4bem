#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 21:26:29 2023

@author: cghiaus
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem


# Physical characteristics
# ========================
concrete = {'Conductivity': 1.400,      # W/(m⋅K)
            'Density': 2300.0,          # kg/m³
            'Specific heat': 880,       # J/(kg⋅K)
            'Width': 0.2}               # m

insulation = {'Conductivity': 0.040,    # W/(m⋅K)
              'Density': 16.0,          # kg/m³
              'Specific heat': 1210,    # J/(kg⋅K)
              'Width': 0.08}            # m

wall = pd.DataFrame.from_dict({'Layer_1': concrete,
                               'Layer_2': insulation},
                              orient='index')

air = {'Density': 1.2,                  # kg/m³
       'Specific heat': 1000}           # J/(kg⋅K)

# convection coefficients, W/(m²·K)
h = pd.DataFrame([{'in': 4., 'out': 10.}], index=['h'])

S_wall = 3 * 3      # m², wall surface area
V_air = 3 * 3 * 3   # m³, indoor air volume

# Resistances
# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_wall)  # K/W
# convection
R_cv = 1 / (h * S_wall)     # K/W

# Capacities
C_wall = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_wall
C_air = air['Density'] * air['Specific heat'] * V_air

# Differential-algebraic equations (DAE)
# =====================================
# number of temperature nodes and flow branches
no_θ = no_q = 7

# Conductance matrix
R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Layer_1'] / 8
R[1] = R[2] = R[3] = R_cd['Layer_1'] / 4
R[4] = R_cd['Layer_1'] / 8 + R_cd['Layer_2'] / 4
R[5] = R_cd['Layer_2'] / 2
R[6] = R_cd['Layer_2'] / 4 + R_cv['in']
G = np.diag(np.reciprocal(R))

# Capacity matrix
C = np.zeros(no_θ)
C[0] = C[1] = C[2] = C[3] = C_wall['Layer_1'] / 4
C[4] = C[5] = C_wall['Layer_2'] / 2
C[6] = C_air
C = np.diag(C)

# Arc-node incidence matrix
A = np.eye(no_q, no_θ + 1)
A = -np.diff(A, n=1, axis=1)

# Input vectors
b = np.zeros(no_q)  # temperatures
f = np.zeros(no_θ)  # flow rates

# Steady-state solution from DAE)
b[0] = 1
θ_steady_To = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
np.set_printoptions(precision=3)
print('When To = 1°C, temperatures in steady-state are:', θ_steady_To, '°C')
print(f'The indoor temperature is: {θ_steady_To[-1]:.3f} °C')

b[0] = 0
f[-1] = 1
θ_steady_Qh = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
print('When Qh = 1W, temperatures in steady-state are:', θ_steady_Qh, '°C')
print(f'The indoor temperature is: {θ_steady_Qh[-1]:.3f} °C')

# State-space representation
# ==========================
# State matrix
As = -np.linalg.inv(C) @ A.T @ G @ A

# Input matrix
Bs = np.linalg.inv(C) @ np.block([A.T @ G, np.eye(no_θ)])
# Select columns for which the input vector is not zero
# 1st for To and last for Qh
Bs = Bs[:, [0, -1]]

# Output matrix
Cs = np.zeros((1, no_θ))
# output: last temperature node
Cs[:, -1] = 1

# Feedthrough (or feedforward) matrix
Ds = np.zeros(Bs.shape[1])

# Eigenvalues
λ = np.linalg.eig(As)[0]    # minimum eigenvalue of matrix A
max_Δt = min(-2 / λ)

np.set_printoptions(precision=1)
print('Time constants: \n', -1 / λ, 's \n')
print('3 x Time constants: \n', -3 / λ, 's \n')
print(f'Max time step Δt = {max_Δt:.2f} s')

# time step
Δt = np.floor(max_Δt / 60) * 60   # s
print(f'Δt = {Δt} s')

# settling time
t_settle = 4 * max(-1 / λ)
print(f'Settling time: {t_settle:.2f} s = {t_settle / 3600:.2f} h')

# Time integration
# ================

# Step input
# ----------

# number of time steps
n = int(np.ceil(t_settle / Δt))
# time vector
t = np.arange(0, n * Δt, Δt)
pd.DataFrame(t, columns=['time'])

# outdoor temperature
# -------------------
u = np.block([[np.ones([1, n])],    # To = [1, 1, ... , 1]
              [np.zeros([1, n])]])  # Qh = [0, 0, ... , 0]

# initial values for temperatures obtained by explicit and implicit Euler
θ_exp = np.zeros([no_θ, t.shape[0]])
θ_imp = np.zeros([no_θ, t.shape[0]])

# time integration: Euler explicit & implicit
for k in range(t.shape[0] - 1):
    θ_exp[:, k + 1] = (np.eye(no_θ) + Δt * As) @\
        θ_exp[:, k] + Δt * Bs @ u[:, k]
    θ_imp[:, k + 1] = np.linalg.inv(np.eye(no_θ) - Δt * As) @\
        (θ_imp[:, k] + Δt * Bs @ u[:, k])

# plot results
fig, ax = plt.subplots()
ax.plot(t / 3600, θ_exp[-1, :], t / 3600, θ_imp[-1, :])
ax.set(xlabel='Time [h]',
       ylabel='Air temperature [°C]',
       title='Step input: $T_o$')
ax.legend(['Explicit', 'Implicit'])
plt.show()

# indoor heat flow rate
# ---------------------
u = np.block([[np.zeros([1, n])],   # To = [0, 0, ... , 0]
              [np.ones([1, n])]])   # Qh = [1, 1, ... , 1]

# initial values for temperatures obtained by explicit and implicit Euler
θ_exp = np.zeros([no_θ, t.shape[0]])
θ_imp = np.zeros([no_θ, t.shape[0]])

# time integration: Euler explicit & implicit
for k in range(t.shape[0] - 1):
    θ_exp[:, k + 1] = (np.eye(no_θ) + Δt * As) @\
        θ_exp[:, k] + Δt * Bs @ u[:, k]
    θ_imp[:, k + 1] = np.linalg.inv(np.eye(no_θ) - Δt * As) @\
        (θ_imp[:, k] + Δt * Bs @ u[:, k])

# plot results
fig, ax = plt.subplots()
ax.plot(t / 3600, θ_exp[-1, :], t / 3600, θ_imp[-1, :])
ax.set(xlabel='Time [h]',
       ylabel='Air temperature [°C]',
       title='Step input: $T_o$')
ax.legend(['Explicit', 'Implicit'])
plt.show()

# Simulation with outdoor temperature from weather data
# -----------------------------------------------------
# Outdoor temperature from weather data
filename = './weather_data/FRA_Lyon.074810_IWEC.epw'
start_date = '2000-04-10'
end_date = '2000-05-15'

[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data
weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather[(
    weather.index >= start_date) & (
    weather.index < end_date)]

# time vector for weather data at 1 h time step
tw = np.arange(0, 3600 * weather.shape[0], 3600)

# time vector for simulation
t = np.arange(0, 3600 * weather.shape[0], Δt)

# resample outdoor temperature at timestep Δt
θ_out = np.interp(t, tw, weather['temp_air'])

# input vector [To, Qh]
u = np.block([[θ_out],
             [np.zeros(θ_out.shape[0])]])

# initial coditions
θ_exp = 20 * np.ones([no_θ, t.shape[0]])
θ_imp = 20 * np.ones([no_θ, t.shape[0]])

# time integration: Euler explicit & implicit
for k in range(u.shape[1] - 1):
    θ_exp[:, k + 1] = (np.eye(no_θ) + Δt * As) @\
        θ_exp[:, k] + Δt * Bs @ u[:, k]
    θ_imp[:, k + 1] = np.linalg.inv(np.eye(no_θ) - Δt * As) @\
        (θ_imp[:, k] + Δt * Bs @ u[:, k])

# plot results
fig, ax = plt.subplots()
ax.plot(t / 3600 / 24, θ_out, label='Outdoor temperature')
ax.plot(t / 3600 / 24, θ_exp[-1, :], label='Indoor temperature')
ax.set(xlabel='Time [days]',
       ylabel='Air temperature [°C]',
       title='Explicit Euler')
ax.legend()
plt.show()

# Simulation with outdoor temperature from weather data with Pandas
# -----------------------------------------------------------------
start_date = '2000-04-10'
end_date = '2000-05-15'

# read data and keep air temperature
filename = './weather_data/FRA_Lyon.074810_IWEC.epw'
[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air"]].copy()
del data

# replace years with year 2000 and select time interval
weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather.loc[start_date:end_date]

# resample weather data
data = weather['temp_air']
data = data.resample(str(Δt) + 'S').interpolate(method='linear')
data = data.rename('To').to_frame()

# indoor auxiliary heat
data['Qa'] = 0 * np.ones(data.shape[0])

# input vector
u = data[['To', 'Qa']]

# initial conditions
θ_exp = 20 * np.ones([As.shape[0], u.shape[0]])
θ_imp = 20 * np.ones([As.shape[0], u.shape[0]])

# time integration: Euler explicit & implicit
n_states = As.shape[0]
I = np.eye(n_states)

for k in range(u.shape[0] - 1):
    θ_exp[:, k + 1] = (I + Δt * As) @ θ_exp[:, k]\
        + Δt * Bs @ u.iloc[k, :]
    θ_imp[:, k + 1] = np.linalg.inv(I - Δt * As) @\
        (θ_imp[:, k] + Δt * Bs @ u.iloc[k, :])

data['θi_exp'] = θ_exp[-1, :]
data['θi_imp'] = θ_imp[-1, :]

ax = data[['To', 'θi_exp']].plot()
# data[['To', 'θi_exp', 'θi_imp']].plot()
ax.legend(['Outdoor temperature', 'Indoor temperature'])
ax.set(xlabel='Time',
       ylabel='Air temperature [°C]',
       title='Explicit Euler')
plt.show()
