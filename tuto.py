#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 19:05:01 2021

@author: cghiaus
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem

# Radiative properties
# --------------------
""" concrete EngToolbox Emissivity Coefficient Materials """
ε_wLW = 0.9     # long wave wall emmisivity
""" grey to dark surface EngToolbox,
    Absorbed Solar Radiation by Surface Color """
α_wSW = 0.2     # absortivity white surface

""" Glass, pyrex EngToolbox Absorbed Solar Radiation bySurface Color """
ε_gLW = 0.9     # long wave glass emmisivity

""" EngToolbox Optical properties of some typical glazing mat
    Window glass """
τ_gSW = 0.83    # short wave glass transmitance

α_gSW = 0.1     # short wave glass absortivity

σ = 5.67e-8     # W/m².K⁴ Stefan-Bolzmann constant
Fwg = 1 / 5     # view factor wall - glass
Tm = 20 + 273   # mean temp for radiative exchange

# Geometry
l = 3                       # m length of the cubic room
Va = l**3                   # m³ volume of air
ACH = 1                     # air changes per hour
Va_dot = ACH * Va / 3600    # m³/s air infiltration

# Thermophyscal properties
# ------------------------
air = {'Density': 1.2,                      # kg/m³
       'Specific heat': 1000}               # J/kg.K

""" Incropera et al. (2011) Fundamantals of heat and mass transfer, 7 ed,
    Table A3,
        concrete (stone mix) p. 993
        insulation polystyrene extruded (R-12) p.990
        glass plate p.993"""
wall = {'Conductivity': [1.4, 0.027, 1.4],      # W/m.K
        'Density': [2300, 55, 2500],            # kg/m³
        'Specific heat': [880, 1210, 750],      # J/kg.K
        'Width': [0.2, 0.08, 0.004],
        'Surface': [5 * l**2, 5 * l**2, l**2],  # m²
        'Meshes': [1, 1, 1]}                      # number of meshes
wall = pd.DataFrame(wall, index=['Concrete', 'Insulation', 'Glass'])

# convection coefficients, W/m² K
h = pd.DataFrame([{'in': 4., 'out': 10}])


def thermal_circuit(Kp):
    """
    Parameters
    ----------
    Kp : float
        Gain of the proportional controller.

    Returns
    -------
    TCa : dict
        {'A': A, 'G': G, 'b': b, 'C': C, 'f': f, 'y': y}.

    DAE Model :
        C * θ_dot = (A.T * G * A) * θ + A.T * G * b + f
        y:  nodes of interest
    """
    # Thermal conductances
    # Conduction
    G_cd = wall['Conductivity'] / wall['Width'] * wall['Surface']

    # Convection
    Gw = h * wall['Surface'][0]     # wall
    Gg = h * wall['Surface'][2]     # glass

    # Long-wave radiation exchnage
    GLW1 = ε_wLW / (1 - ε_wLW) * wall['Surface']['Insulation'] * 4 * σ * Tm**3
    GLW2 = Fwg * wall['Surface']['Insulation'] * 4 * σ * Tm**3
    GLW3 = ε_gLW / (1 - ε_gLW) * wall['Surface']['Glass'] * 4 * σ * Tm**3
    # long-wave exg. wall-glass
    GLW = 1 / (1 / GLW1 + 1 / GLW2 + 1 / GLW3)

    # ventilation & advection
    Gv = Va_dot * air['Density'] * air['Specific heat']

    # glass: convection outdoor & conduction
    Ggs = float(1 / (1 / Gg['out'] + 1 / (2 * G_cd['Glass'])))

    # Thermal capacities
    Capacity = wall['Density'] * wall['Specific heat'] *\
        wall['Surface'] * wall['Width']
    Capacity['Air'] = air['Density'] * air['Specific heat'] * Va

    # Thermal network
    # ---------------
    # Dissembled circuits
    # TCd0:  Concrete and insulation wall  (in red)
    nq = 1 + 2 * (wall['Meshes']['Concrete'] + wall['Meshes']['Insulation'])
    nt = 1 + 2 * (wall['Meshes']['Concrete'] + wall['Meshes']['Insulation'])

    A = np.eye(nq + 1, nt)
    A = -np.diff(A, 1, 0).T

    nc = wall['Meshes']['Concrete']
    ni = wall['Meshes']['Insulation']
    Gcm = 2 * nc * [G_cd['Concrete']]
    Gcm = 2 * nc * np.array(Gcm)
    Gim = 2 * ni * [G_cd['Insulation']]
    Gim = 2 * ni * np.array(Gim)
    G = np.diag(np.hstack([Gw['out'], Gcm, Gim]))

    b = np.zeros(nq)
    b[0] = 1

    Ccm = Capacity['Concrete'] / nc * np.mod(range(0, 2 * nc), 2)
    Cim = Capacity['Insulation'] / ni * np.mod(range(0, 2 * ni), 2)
    C = np.diag(np.hstack([Ccm, Cim, 0]))

    f = np.zeros(nt)
    f[0] = f[-1] = 1

    y = np.zeros(nt)

    TCd0 = {'A': A, 'G': G, 'b': b, 'C': C, 'f': f, 'y': y}

    # TCd1: Indoor air (in blue)
    A = np.array([[-1, 1, 0],
                  [-1, 0, 1],
                  [0, -1, 1]])
    G = np.diag(np.hstack([GLW, Gw['in'], h['in'] * wall['Surface']['Glass']]))
    b = np.zeros(3)
    C = np.diag([0, 0, Capacity['Air'] / 2])
    f = np.array([1, 0, 1])
    y = np.array([0, 0, 1])
    TCd1 = {'A': A, 'G': G, 'b': b, 'C': C, 'f': f, 'y': y}

    # TCd2: Glass (in green)
    A = np.array([[1, 0],
                  [-1, 1]])
    Ggo = h['out'] * wall['Surface']['Glass']
    Ggs = 1 / (1 / Ggo + 1 / (2 * G_cd['Glass']))
    G = np.diag(np.hstack([Ggs, 2 * G_cd['Glass']]))
    b = np.array([1, 0])
    C = np.diag([Capacity['Glass'], 0])
    f = np.array([1, 0])
    y = np.array([0, 0])
    TCd2 = {'A': A, 'G': G, 'b': b, 'C': C, 'f': f, 'y': y}

    # TCd3: air infiltration and controller (in purple)
    A = np.array([[1],
                  [1]])
    G = np.diag(np.hstack([Gv, Kp]))
    b = np.array([1, 1])
    C = np.array([Capacity['Air'] / 2])
    f = 1
    y = 1
    TCd3 = {'A': A, 'G': G, 'b': b, 'C': C, 'f': f, 'y': y}

    TCd = {'0': TCd0,
           '1': TCd1,
           '2': TCd2,
           '3': TCd3}

    AssX = np.array([[0, nt - 1, 1, 0],
                     [1, 1, 2, 1],
                     [1, 2, 3, 0]])

    TCa = dm4bem.TCAss(TCd, AssX)
    return TCa


def step_response(duration, dt, As, Bs, Cs, Ds):
    # Step response
    # -------------
    # number of steps
    n = int(np.floor(duration / dt))

    t = np.arange(0, n * dt, dt)    # time

    # Vectors of state and input (in time)
    n_tC = As.shape[0]              # no of state variables
    # u = [To To To Tsp Phio Phii Qaux Phia]
    u = np.zeros([8, n])
    u[0:3, :] = np.ones([3, n])
    # initial values for temperatures obtained by explicit and implicit Euler
    temp_exp = np.zeros([n_tC, t.shape[0]])
    temp_imp = np.zeros([n_tC, t.shape[0]])

    I = np.eye(n_tC)
    for k in range(n - 1):
        temp_exp[:, k + 1] = (I + dt * As) @\
            temp_exp[:, k] + dt * Bs @ u[:, k]
        temp_imp[:, k + 1] = np.linalg.inv(I - dt * As) @\
            (temp_imp[:, k] + dt * Bs @ u[:, k])

    y_exp = Cs @ temp_exp + Ds @  u
    y_imp = Cs @ temp_imp + Ds @  u

    fig, ax = plt.subplots()
    ax.plot(t / 3600, y_exp.T, label='Euler explicit')
    ax.plot(t / 3600, y_imp.T, label='Euler implicit')
    ax.set(xlabel='Time [h]',
           ylabel='$T_i$ [°C]',
           title='Step input: To = 1°C')
    ax.legend(loc='lower right')


def P_control(filename, start_date, end_date, dt,
              As, Bs, Cs, Ds, Kp):
    # Read weather data from Energyplus .epw file
    [data, meta] = dm4bem.read_epw(filename, coerce_year=None)
    weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
    del data
    weather.index = weather.index.map(lambda t: t.replace(year=2000))
    weather = weather[(weather.index >= start_date) & (
        weather.index < end_date)]

    # Solar radiation on a tilted surface
    surface_orientation = {'slope': 90,
                           'azimuth': 0,
                           'latitude': 45}
    albedo = 0.2
    rad_surf1 = dm4bem.sol_rad_tilt_surf(weather, surface_orientation, albedo)
    rad_surf1['Φt1'] = rad_surf1.sum(axis=1)

    # Interpolate weather data for time step dt
    data = pd.concat([weather['temp_air'], rad_surf1['Φt1']], axis=1)
    data = data.resample(str(dt) + 'S').interpolate(method='linear')
    data = data.rename(columns={'temp_air': 'To'})

    # Indoor temperature set-point
    data['Ti'] = 20 * np.ones(data.shape[0])

    # Indoor auxiliary heat flow rate
    data['Qa'] = 0 * np.ones(data.shape[0])

    # time
    t = dt * np.arange(data.shape[0])

    u = pd.concat([data['To'], data['To'], data['To'], data['Ti'],
                   α_wSW * wall['Surface']['Concrete'] * data['Φt1'],
                   τ_gSW * α_wSW * wall['Surface']['Glass'] * data['Φt1'],
                   data['Qa'],
                   α_gSW * wall['Surface']['Glass'] * data['Φt1']], axis=1)

    # initial values for temperatures
    temp_exp = 20 * np.ones([As.shape[0], u.shape[0]])

    # integration in time
    I = np.eye(As.shape[0])
    for k in range(u.shape[0] - 1):
        temp_exp[:, k + 1] = (I + dt * As) @ temp_exp[:, k]\
            + dt * Bs @ u.iloc[k, :]
    # Indoor temperature
    y_exp = Cs @ temp_exp + Ds @ u.to_numpy().T
    # HVAC heat flow
    q_HVAC = Kp * (data['Ti'] - y_exp[0, :])

    # plot indoor and outdoor temperature
    fig, axs = plt.subplots(2, 1)
    axs[0].plot(t / 3600, y_exp[0, :], label='$T_{indoor}$')
    axs[0].plot(t / 3600, data['To'], label='$T_{outdoor}$')
    axs[0].set(xlabel='Time [h]',
               ylabel='Temperatures [°C]',
               title='Simulation for weather')
    axs[0].legend(loc='upper right')

    # plot total solar radiation and HVAC heat flow
    axs[1].plot(t / 3600, q_HVAC, label='$q_{HVAC}$')
    axs[1].plot(t / 3600, data['Φt1'], label='$Φ_{total}$')
    axs[1].set(xlabel='Time [h]',
               ylabel='Heat flows [W]')
    axs[1].legend(loc='upper right')

    fig.tight_layout()

    return y_exp, q_HVAC


def switch_models(filename, start_date, end_date, dt,
                 Af, Bf, Cf, Df, Kpf,
                 Ac, Bc, Cc, Dc, Kpc,
                 Tisp, DeltaT):
    """
    Use of two models, one in free-running and one perfect controller

    Parameters
    ----------
    filename : TYPE
        Weather file.
    start_date : TYPE
        DESCRIPTION.
    end_date : TYPE
        DESCRIPTION.
    dt : TYPE
        Integration time step.
    Af, Bf, Cf, Df : TYPE
        State-space model for free-running.
    Kpf : TYPE
        Controller gain for free-running (Kpf -> 0).
    Ac, Bc, Cc, Dc : TYPE
        State-space model for perfect control..
    Kpc : TYPE
        Controller gain for perfect P-control (Kpc -> infinity).
    Tisp : TYPE
        Indoor temperature set point.
    DeltaT : TYPE
        Dead-band (accepted switch) temperature.

    Returns
    -------
    None.

    """
    # Read weather data from Energyplus .epw file
    [data, meta] = dm4bem.read_epw(filename, coerce_year=None)
    weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
    del data
    weather.index = weather.index.map(lambda t: t.replace(year=2000))
    weather = weather[(weather.index >= start_date) & (
        weather.index < end_date)]

    # Solar radiation on a tilted surface
    surface_orientation = {'slope': 90,
                           'azimuth': 0,
                           'latitude': 45}
    albedo = 0.2
    rad_surf1 = dm4bem.sol_rad_tilt_surf(weather, surface_orientation, albedo)
    rad_surf1['Φt1'] = rad_surf1.sum(axis=1)

    # Interpolate weather data for time step dt
    data = pd.concat([weather['temp_air'], rad_surf1['Φt1']], axis=1)
    data = data.resample(str(dt) + 'S').interpolate(method='linear')
    data = data.rename(columns={'temp_air': 'To'})

    # Indoor temperature set-point
    data['Ti'] = 20 * np.ones(data.shape[0])

    # Indoor auxiliary heat flow rate
    data['Qa'] = 0 * np.ones(data.shape[0])

    # time
    t = dt * np.arange(data.shape[0])

    u = pd.concat([data['To'], data['To'], data['To'], data['Ti'],
                   α_wSW * wall['Surface']['Concrete'] * data['Φt1'],
                   τ_gSW * α_wSW * wall['Surface']['Glass'] * data['Φt1'],
                   data['Qa'],
                   α_gSW * wall['Surface']['Glass'] * data['Φt1']], axis=1)

    # initial values for temperatures
    temp_exp = 0 * np.ones([Af.shape[0], u.shape[0]])
    Tisp = Tisp * np.ones(u.shape[0])
    y = np.zeros(u.shape[0])
    y[0] = Tisp[0]
    q_HVAC = 0 * np.ones(u.shape[0])

    # integration in time
    I = np.eye(Af.shape[0])
    for k in range(u.shape[0] - 1):
        if y[k] < Tisp[k] or y[k] > DeltaT + Tisp[k]:
            temp_exp[:, k + 1] = (I + dt * Ac) @ temp_exp[:, k]\
                + dt * Bc @ u.iloc[k, :]
            y[k + 1] = Cc @ temp_exp[:, k + 1] + Dc @ u.iloc[k + 1]
            q_HVAC[k + 1] = Kpc * (Tisp[k + 1] - y[k + 1])
        else:
            temp_exp[:, k + 1] = (I + dt * Af) @ temp_exp[:, k]\
                + dt * Bf @ u.iloc[k, :]
            y[k + 1] = Cf @ temp_exp[:, k + 1] + Df @ u.iloc[k]
            q_HVAC[k + 1] = 0

    # plot indoor and outdoor temperature
    fig, axs = plt.subplots(2, 1)
    axs[0].plot(t / 3600, y, label='$T_{indoor}$')
    axs[0].plot(t / 3600, data['To'], label='$T_{outdoor}$')
    axs[0].set(xlabel='Time [h]',
               ylabel='Temperatures [°C]',
               title='Simulation for weather')
    axs[0].legend(loc='upper right')

    # plot total solar radiation and HVAC heat flow
    axs[1].plot(t / 3600, q_HVAC, label='$q_{HVAC}$')
    axs[1].plot(t / 3600, data['Φt1'], label='$Φ_{total}$')
    axs[1].set(xlabel='Time [h]',
               ylabel='Heat flows [W]')
    axs[1].legend(loc='upper right')
    plt.ylim(-1500, 3000)
    fig.tight_layout()

    return y, q_HVAC


def heat(filename, start_date, end_date, dt,
         Af, Bf, Cf, Df, Kpf,
         Ac, Bc, Cc, Dc, Kpc,
         Tisp, DeltaT):
    """
    Use of two models, one in free-running and one perfect controller

    Parameters
    ----------
    filename : TYPE
        Weather file.
    start_date : TYPE
        DESCRIPTION.
    end_date : TYPE
        DESCRIPTION.
    dt : TYPE
        Integration time step.
    Af, Bf, Cf, Df : TYPE
        State-space model for free-running.
    Kpf : TYPE
        Controller gain for free-running (Kpf -> 0).
    Ac, Bc, Cc, Dc : TYPE
        State-space model for perfect control..
    Kpc : TYPE
        Controller gain for perfect P-control (Kpc -> infinity).
    Tisp : TYPE
        Indoor temperature set point.
    DeltaT : TYPE
        Dead-band (accepted switch) temperature.

    Returns
    -------
    None.

    """
    # Read weather data from Energyplus .epw file
    [data, meta] = dm4bem.read_epw(filename, coerce_year=None)
    weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
    del data
    weather.index = weather.index.map(lambda t: t.replace(year=2000))
    weather = weather[(weather.index >= start_date) & (
        weather.index < end_date)]

    # Solar radiation on a tilted surface
    surface_orientation = {'slope': 90,
                           'azimuth': 0,
                           'latitude': 45}
    albedo = 0.2
    rad_surf1 = dm4bem.sol_rad_tilt_surf(weather, surface_orientation, albedo)
    rad_surf1['Φt1'] = rad_surf1.sum(axis=1)

    # Interpolate weather data for time step dt
    data = pd.concat([weather['temp_air'], rad_surf1['Φt1']], axis=1)
    data = data.resample(str(dt) + 'S').interpolate(method='linear')
    data = data.rename(columns={'temp_air': 'To'})

    # Indoor temperature set-point
    data['Ti'] = 20 * np.ones(data.shape[0])

    # Indoor auxiliary heat flow rate
    data['Qa'] = 0 * np.ones(data.shape[0])

    # time
    t = dt * np.arange(data.shape[0])

    u = pd.concat([data['To'], data['To'], data['To'], data['Ti'],
                   α_wSW * wall['Surface']['Concrete'] * data['Φt1'],
                   τ_gSW * α_wSW * wall['Surface']['Glass'] * data['Φt1'],
                   data['Qa'],
                   α_gSW * wall['Surface']['Glass'] * data['Φt1']], axis=1)

    # initial values
    temp_exp = 0 * np.ones([Af.shape[0], u.shape[0]])
    temp_exp[:, 0] = data['Ti'][0]

    y = np.zeros(u.shape[0])
    y[0] = data['Ti'][0]
    q_HVAC = 0 * np.ones(u.shape[0])

    # integration in time
    I = np.eye(Af.shape[0])

    for k in range(u.shape[0] - 1):
        if y[k] < data['Ti'][k]:
            temp_exp[:, k + 1] = (I + dt * Ac) @ temp_exp[:, k]\
                + dt * Bc @ u.iloc[k, :]
            y[k + 1] = Cc @ temp_exp[:, k + 1] + Dc @ u.iloc[k]
            q_HVAC[k + 1] = Kpc * (data['Ti'][k + 1] - y[k + 1])
        else:
            temp_exp[:, k + 1] = (I + dt * Af) @ temp_exp[:, k]\
                + dt * Bf @ u.iloc[k, :]
            y[k + 1] = Cf @ temp_exp[:, k + 1] + Df @ u.iloc[k]
            q_HVAC[k + 1] = 0

    # plot indoor and outdoor temperature
    fig, axs = plt.subplots(2, 1)
    axs[0].plot(t / 3600, y, label='$T_{indoor}$')
    axs[0].plot(t / 3600, data['To'], label='$T_{outdoor}$')
    axs[0].set(xlabel='Time [h]',
               ylabel='Temperatures [°C]',
               title='Simulation for weather')
    axs[0].legend(loc='upper right')

    # plot total solar radiation and HVAC heat flow
    axs[1].plot(t / 3600, q_HVAC, label='$q_{HVAC}$')
    axs[1].plot(t / 3600, data['Φt1'], label='$Φ_{total}$')
    axs[1].set(xlabel='Time [h]',
               ylabel='Heat flows [W]')
    axs[1].legend(loc='upper right')
    plt.ylim(-1500, 3000)
    fig.tight_layout()

    return y, q_HVAC


def inputs(filename, start_date, end_date, dt,
           As, Bs, Cs, Ds, Kp, Tisp):
    # Read weather data from Energyplus .epw file
    [data, meta] = dm4bem.read_epw(filename, coerce_year=None)
    weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
    del data
    weather.index = weather.index.map(lambda t: t.replace(year=2000))
    weather = weather[(weather.index >= start_date) & (
        weather.index < end_date)]

    # Solar radiation on a tilted surface
    surface_orientation = {'slope': 90,
                           'azimuth': 0,
                           'latitude': 45}
    albedo = 0.2
    rad_surf1 = dm4bem.sol_rad_tilt_surf(weather, surface_orientation, albedo)
    rad_surf1['Φt1'] = rad_surf1.sum(axis=1)

    # Interpolate weather data for time step dt
    data = pd.concat([weather['temp_air'], rad_surf1['Φt1']], axis=1)
    data = data.resample(str(dt) + 'S').interpolate(method='linear')
    data = data.rename(columns={'temp_air': 'To'})

    # Indoor temperature set-point
    data['Ti'] = Tisp * np.ones(data.shape[0])

    # Indoor auxiliary heat flow rate
    data['Qa'] = 0 * np.ones(data.shape[0])

    # time
    t = dt * np.arange(data.shape[0])

    u = pd.concat([data['To'], data['To'], data['To'], data['Ti'],
                   α_wSW * wall['Surface']['Concrete'] * data['Φt1'],
                   τ_gSW * α_wSW * wall['Surface']['Glass'] * data['Φt1'],
                   data['Qa'],
                   α_gSW * wall['Surface']['Glass'] * data['Φt1']], axis=1)

    return t, u, data


def plot_results(t, y, temp_exp, q_HVAC, data):
    """

    Parameters
    ----------
    t : Array of floats
        Time in seconds.
    y : Array of floats
        Indoor temperature, °C.
    temp_exp : TYPE
        DESCRIPTION.
    q_HVAC : Array of floats
        Heat flow rate of HVAC system, W.
    data : DataFrame
        Index: time at dt
        To:     outdoor temperature, °C
        Φt1:    total solar radiation, W / m2;
        Ti:     indoor setpoint temperature, °C;
        Qa:     auxiliary heat flow rate, W.

    Returns
    -------
    None.

    """
    fig, axs = plt.subplots(2, 1)
    axs[0].plot(t / 3600, y, label='$T_{indoor}$')
    axs[0].plot(t / 3600, data['To'], label='$T_{outdoor}$')
    # axs[0].plot(t / 3600, temp_exp[2], label='$T_{ss}$')
    axs[0].set(xlabel='Time [h]',
               ylabel='Temperatures [°C]',
               title='Simulation for weather')
    axs[0].legend(loc='upper right')
    axs[0].set_ylim([-2, 40])

    # plot total solar radiation and HVAC heat flow
    axs[1].plot(t / 3600, q_HVAC, label='$q_{HVAC}$')
    axs[1].plot(t / 3600, data['Φt1'], label='$Φ_{total}$')
    axs[1].set(xlabel='Time [h]',
               ylabel='Heat flows [W]')
    axs[1].legend(loc='upper right')
    axs[1].set_ylim([-1000, 2000])
    fig.tight_layout()
    plt.draw()
    plt.show()
