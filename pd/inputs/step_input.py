#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:32:24 2023

@author: cghiaus

Simulate the model from /dm4bem/03CubicBuilding.ipynb for step inputs

Steps:
    - from /bldg -> disassambled TCd -> assembled TC -> state-space
    - create step inputs: To, Ti_sp, Φo, Φi, Qa, Φa
    - obtain input vector u = [TO, To, To, Ti_sp, Φo, Φi, Qa, Φa]
    - simulate ss-model for u
"""
import numpy as np
import pandas as pd
import pd_dm4bem

# Obtain state-space representation
# =================================
# Disassembled thermal circuits
folder_path = "bldg"
TCd = pd_dm4bem.bldg2TCd(folder_path,
                         TC_auto_number=True)

# Assembled thermal circuit using assembly_lists.csv'
ass_lists = pd.read_csv(folder_path + '/assembly_lists.csv')
ass_matrix = pd_dm4bem.assemble_lists2matrix(ass_lists)
TC = pd_dm4bem.assemble_TCd_matrix(TCd, ass_matrix)

# State-space
[As, Bs, Cs, Ds, us] = pd_dm4bem.tc2ss(TC)

# Obtain input vector in time
# ===========================
# Eigenvalue analysis
λ = np.linalg.eig(As)[0]    # eigenvalues of matrix As
λ = np.sort(λ)

print('Time constants:')
print([f'{T:.2f} s' for T in -1 / λ])

dt_max = 2 * min(-1. / λ)
print(f'\nMaximum time step: {dt_max:.2f} s = {dt_max / 60:.2f} min')

t_settle = 4 * max(-1. / λ)
print(f'Minimum settling time: \
{t_settle:.0f} s = \
{t_settle / 60:.1f} min = \
{t_settle / 3600:.2f} h = \
{t_settle / (3600 * 24):.2f} days')

# time step
dt = np.floor(dt_max / 60) * 60   # s
print(f'dt = {dt} s = {dt / 60:.0f} min')

# duration: next multiple of 3600 s that is larger than t_settle
duration = np.ceil(t_settle / 3600) * 3600
print(f'Duration = {duration} s')

# Create input_data_set
# ---------------------
# time vector
n = int(np.floor(duration / dt))    # number of time steps

# Create a DateTimeIndex starting at "00:00:00" with a time step of dt
time = pd.date_range(start="2000-01-01 00:00:00",
                           periods=n, freq=f"{int(dt)}S")

To = 10 * np.ones(n)
Ti_sp = 20 * np.ones(n)
Φa = 0 * np.ones(n)
Qa = Φo = Φi = Φa

data = {'To': To, 'Ti_sp': Ti_sp, 'Qa': Qa, 'Φo': Φo, 'Φi': Φi, 'Φa': Φa}
input_data_set = pd.DataFrame(data, index=time)

# Get input from input_data_set
u = pd_dm4bem.inputs_in_time(us, input_data_set)

# Simulation
# ==========
# Initial conditions
θ0 = 0                      # initial temperatures
θ_exp = pd.DataFrame(index=u.index)
θ_exp[As.columns] = θ0      # Fill θ with initial valeus θ0
θ_imp = θ_exp

I = np.eye(As.shape[0])     # identity matrix

for k in range(n - 1):
    θ_exp.iloc[k + 1] = (I + dt * As)\
        @ θ_exp.iloc[k] + dt * Bs @ u.iloc[k]
    θ_imp.iloc[k + 1] = np.linalg.inv(I - dt * As)\
        @ (θ_imp.iloc[k] + dt * Bs @ u.iloc[k])

# outputs
y_exp = (Cs @ θ_exp.T + Ds @  u.T).T
y_imp = (Cs @ θ_imp.T + Ds @  u.T).T

# Plot resukts
y = pd.concat([y_exp, y_imp], axis=1, keys=['Explicit', 'Implicit'])
# Flatten the two-level column labels into a single level
y.columns = y.columns.get_level_values(0)
y.plot()
