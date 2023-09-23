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
import matplotlib.pyplot as plt
import pd_dm4bem

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

# time vector
n = int(np.floor(duration / dt))    # number of time steps

# Create a DateTimeIndex starting at "00:00:00" with a time step of dt
time = pd.date_range(start="2000-01-01 00:00:00",
                           periods=n, freq=f"{int(dt)}S")

# # Create the Pandas Series 'To' with constant values
# To = pd.Series(10 * np.ones(n), index=time)
# Ti_sp = pd.Series(20 * np.ones(n), index=time)
# Φa = pd.Series(0 * np.ones(n), index=time)
# Qa = Φo = Φi = Φa


# # Create the DataFrame 'u' based on the 'us' Series
# u = pd.DataFrame({col: globals()[us[col]] for col in us.index})

# # is achieving
# # u = pd.DataFrame({'c1_q0': To, 'c2_q0': To, 'c3_q0': Ti_sp, 'ow0_q0': To,
# #                   'c1_θ0': Φa, 'c2_θ0': Qa, 'ow0_θ0': Φo, 'ow0_θ4': Φi})
# # from
# # us =
# # c1_q0        To
# # c2_q0        To
# # c3_q0     Ti_sp
# # ow0_q0       To
# # c1_θ0        Φa
# # c2_θ0        Qa
# # ow0_θ0       Φo
# # ow0_θ4       Φi


To = 10 * np.ones(n)
Ti_sp = 20 * np.ones(n)
Φa = 0 * np.ones(n)
Qa = Φo = Φi = Φa

data = {'To': To, 'Ti_sp': Ti_sp, 'Qa': Qa, 'Φo': Φo, 'Φi': Φi, 'Φa': Φa}
u_set = pd.DataFrame(data, index=time)

# # Construct the DataFrame
# u1 = pd.DataFrame({'c1_q0': u_set['To'], 'c2_q0': u_set['To'],
#                    'c3_q0': u_set['Ti_sp'], 'ow0_q0': u_set['To'],
#                    'c1_θ0': u_set['Φa'], 'c2_θ0': u_set['Qa'],
#                    'ow0_θ0': u_set['Φo'], 'ow0_θ4': u_set['Φi']})

# u1 = pd.DataFrame({'c1_q0': u_set['To'], 'c2_q0': u_set['To']})
# u2 = pd.DataFrame({'c1_q0': u_set[us['c1_q0']], 'c2_q0': u_set['To']})

# u3 = pd.DataFrame({col: u_set[us[col]] for col in us.index})

u4 = pd_dm4bem.inputs_in_time(us, u_set)


# u_list = [u_set[column] for column in u_set.columns]


# # Time integration
# # initial conditions
# n_s = As.shape[0]                      # number of state variables
# θ_exp = np.zeros([n_s, time.shape[0]])    # explicit Euler in time t
# θ_imp = np.zeros([n_s, time.shape[0]])    # implicit Euler in time t

# # time integration
# I = np.eye(n_s)                        # identity matrix

# for k in range(n - 1):
#     θ_exp[:, k + 1] = (I + dt * As) @ θ_exp[:, k]\
#         + dt * Bs @ u.iloc[k, :]
#     θ_imp[:, k + 1] = np.linalg.inv(I - dt * As) @ (θ_imp[:, k]\
#                                                     + dt * Bs @ u.iloc[k, :])
# # outputs
# y_exp = Cs @ θ_exp + Ds @  u.T
# y_imp = Cs @ θ_imp + Ds @  u.T

# fig, ax = plt.subplots()
# ax.plot(time / 3600, y_exp.T, time / 3600, y_imp.T)
# ax.set(xlabel='Time, $t$ / h',
#        ylabel='Temperatue, $θ_i$ / °C',
#        title='Step input: outdoor temperature $T_o$')
# ax.legend(['Explicit', 'Implicit'])
# ax.grid()
# plt.show()
