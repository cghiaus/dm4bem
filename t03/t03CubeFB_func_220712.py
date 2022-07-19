#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 13:47:15 2022
Modified   Tue Jul 12 15:11:26 2022

@author: cghiaus


To refactor:
    - read data from .csv files (material properties)

t03CubeFB Tutorial 3: Cubic buiding with feed-back P-controller

Functional implementation of t03CubeFB.ipynb
Uses pd.DataFrame (with index: time) for results of simulations


Structure diagrame (program architecture):

Program structure:
    - conductance_capacity()
    - <step_plot():
        - thermal_circuit()
        - dm4bem.tc2ss()
        - <step_response():
            - time_step()
            - time_response()
            - <integration_Euler()

    - conductance_capacity():
    - <weather_plot():
        - thermal_circuit()
        - dm4bem.tc2ss()
        - <weather_response():
            - time_step()
            - dm4bem.read_epw()
            - dm4bem.sol_rad_tilt_surf()
            - <integration_Euler()

Bibliography:
    Help Docstring Conventions https://peps.python.org/pep-0257/
    https://en.m.wikipedia.org/wiki/Unified_Modeling_Language

Notes:
    iPython console:
    %matplotlib qt      --> plot in separate window
    %matplotlib inline  --> plot in Plots pan. Options Plots pan: Mute inline
"""
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import sys
import dm4bem


def conductance_capacity(ε_wLW=0.85, α_wSW=0.25,
                         ε_gLW=0.90, τ_gSW=0.30, α_gSW=0.38,
                         Fwg=1 / 5,
                         Tm=20 + 273,
                         l=3,
                         ACH=1,
                         air={'ρ': 1.2,
                              'c': 1000},
                         concrete={'λ': 1.4,
                                   'ρ': 2300,
                                   'c': 880,
                                   'w': 0.2,
                                   'n': 1},
                         insulation={'λ': 0.027,
                                     'ρ': 55,
                                     'c': 1210,
                                     'w': 0.08,
                                     'n': 1},
                         glass={'λ': 1.4,
                                'ρ': 2500,
                                'c': 750,
                                'w': 0.004,
                                'n': 1},
                         h={'in': 8,
                            'out': 25},
                         Kp=1e-3):
    """
    Conductances and capacities as a function of thermophysical
    & radiative properties and geometrical dimensions for a
    cubic building having 5 identical walls (2 materials) and 1 glass wall.

    Parameters
    ----------
    Radiative properties (Engineering Toolbox):

    ε_wLW : float, adimensional [0, 1]
        long wave wall emmisivity, 0.85 - grey to dark surface.
    α_wSW : float, adimensional [0, 1]
        short wave absortivity, 0.2 - white surface.
    ε_gLW : float, adimensional [0, 1]
        long wave glass (pyrex) emmisivity.
    ε_wLW : float, adimensional [0, 1]
        long wave wall emmisivity, grey to dark surface.
    α_wSW : float, adimensional [0, 1]
        short wave absortivity, white surface.
    ε_gLW : float, adimensional [0, 1]
        long wave glass (pyrex) emmisivity.

    Optical properties of typical glazing materials:

    τ_gSW : float, adimensional [0, 1]
        short wave glass transmitance of window glass.
    α_gSW : float, adimensional [0, 1]
        short wave glass absortivity of window glass.

    Fwg : float, adimensional [0, 1]
        view factor wall - glass.
    Tm : float, adimensional, K
        mean temperature of surfaces for linearization of radiative exchange.

    Building geometrical dimensions:

    l : float, m
        length of the cubic room.

    Infiltration rate:

    ACH : float, 1/hour
        air changes per hour.

    Thermopysical properties (in SI units):
        - air, dict :
            - ρ : density, kg/m³
            - c : specific heat, J/kg·K
        - concrete, dict :
            - λ : thermal conductivity,  W/m·K
            - ρ : density, kg/m³
            - c : specific heat, J/kg·K
            - w : width of the layer, m
            - n : number of meshes, -,
        - insulation, dict :
            - λ : thermal conductivity,  W/m·K
            - ρ : density, kg/m³
            - c : specific heat, J/kg·K
            - w : width of the layer, m
            - n : number of meshes, -
        - glass, dict :
            - λ : thermal conductivity,  W/m·K
            - ρ : density, kg/m³
            - c : specific heat, J/kg·K
            - w : width of the layer, m
            - n : number of meshes, -

    h : dict, heat convection coeficientsn, W/m²·K :
        - in : float positive, inside surface
        - out : float positive, outside surface

    Kp : float, positif.
        Gain of P-controller :
            - Kp --> 0    free running
            - Kp --> ∞    perfect controller


    Returns
    -------
    G : pandas.Series
    Thermal conductances, W/K :
     - ['cd']:
         - ['concrete']
         - ['insulation']
         - ['glass']
     - ['cv wall']:
         - ['in']
         - ['out'];
     - ['cv glass']:
         - ['in']
         - ['out'];
     - ['LW']
     - ['adv']
     - ['cdv glass']


    C : pandas.Series
    Thermal capacities, J/kg :
        - ['concrete']
        - ['insulation']
        - ['air']

    wall['n'] : int
        Number of meshes in each layer from type(wall) = pandas.DataFrame


    Absortivity × Radiative surface (dict)

    αS: pd.Series
        Abdortivity * Surface for :
            - ['wall_out']  for [α_wSW * wall['Surface']['Concrete']
            - ['wall in']   for [τ_gSW * α_wSW * wall['Surface']['Glass']]
            - ['glass']     for [α_gSW * wall['Surface']['Glass']
    """

    σ = 5.67e-8     # W/(m²·K⁴) Stefan-Bolzmann constant

    wall = {'λ': [concrete['λ'], insulation['λ'], glass['λ']],
            'ρ': [concrete['ρ'], insulation['ρ'], glass['ρ']],
            'c': [concrete['c'], insulation['c'], glass['c']],
            'w': [concrete['w'], insulation['w'], glass['w']],
            'S': [5 * l**2, 5 * l**2, l**2],
            'n': [concrete['n'], insulation['n'], glass['n']]}
    wall = pd.DataFrame(wall, index=['concrete', 'insulation', 'glass'])

    h = pd.Series(h)

    Va = l**3                   # m³ volume of air
    Va_dot = ACH * Va / 3600    # m³/s air infiltration rate

    # Thermal conductances
    #
    # Conduction
    Gw_cd = wall['λ'] / wall['w'] * wall['S']

    # Convection
    Gw_cv = h * wall['S'][0]     # wall
    Gg_cv = h * wall['S'][2]     # glass

    # Long-wave radiation exchange liniearized (Fig. 4)
    GLW1 = ε_wLW / (1 - ε_wLW) * wall['S']['insulation'] * 4 * σ * Tm**3
    GLW2 = Fwg * wall['S']['insulation'] * 4 * σ * Tm**3
    GLW3 = ε_gLW / (1 - ε_gLW) * wall['S']['glass'] * 4 * σ * Tm**3
    # equivalent conductance for long-wave exchange wall-glass
    GLW = 1 / (1 / GLW1 + 1 / GLW2 + 1 / GLW3)

    # Ventilation & advection
    Gv = Va_dot * air['ρ'] * air['c']

    # Glass: convection outdoor & conduction
    Gg_cdv = float(1 / (1 / Gg_cv['out'] + 1 / (2 * Gw_cd['glass'])))

    # Returns
    # -------
    # Thermal conductances
    G = {'cd': Gw_cd,           # G['cd']['concrete'];['insulation'];['glass']
         'cv wall': Gw_cv,      # G['cv wall']['in'];['out']
         'cv glass': Gg_cv,     # G['cv glass']['in'];['out']
         'LW': GLW,             # G['LW']
         'adv': Gv,             # G['adv']
         'cdv glass': Gg_cdv,   # G['cdv glass']
         'Kp': Kp}              # G['Kp']
    G = pd.Series(G)

    # Thermal capacities        # C['concrete'];['insulation'];['air']
    C = wall['ρ'] * wall['c']
    C *= wall['S'] * wall['w']
    C['air'] = air['ρ'] * air['c'] * Va

    # Absortivity × surface of radiation, m²
    αS = {
        'wall out': α_wSW * wall['S']['concrete'],
        'wall in': α_wSW * τ_gSW * wall['S']['glass'],
        'glass': α_gSW * wall['S']['glass']}
    αS = pd.Series(αS)

    return G, C, wall['n'], αS


def thermal_circuit(G, C):
    """
    Thermal circuit (A, G, b, C, f, y) for cubic building (Fig. 3)

    .. math::
        C\\dot{θ} = -A^T G A θ + A^T G b + f

    from conductances & capacities.

    Parameters
    ----------
    G : pandas.Series
        Conductances given by function conductance_capacity()
    C : pandas.Series
        Capacities given by function conductance_capacity()

    Returns
    -------
    A, G, b, C, f, y
        Thermal circuit matrices and vectors numpy.arrays

    """
    nq = 12                             # n° of flow branches
    nθ = 8                              # n° of temperature nodes

    # A: incidence matrix
    A = np.zeros([nq, nθ])              # (n° of branches) × (n° of nodes)
    A[0, 0] = 1                         # branch 0: -> node 0
    A[1, 0], A[1, 1] = -1, 1            # branch 1: node 0 -> node 1
    A[2, 1], A[2, 2] = -1, 1            # branch 2: node 1 -> node 2
    A[3, 2], A[3, 3] = -1, 1            # branch 3: node 2 -> node 3
    A[4, 3], A[4, 4] = -1, 1            # branch 4: node 3 -> node 4
    A[5, 4], A[5, 5] = -1, 1            # branch 5: node 4 -> node 5
    A[6, 4], A[6, 6] = -1, 1            # branch 6: node 4 -> node 6
    A[7, 5], A[7, 6] = -1, 1            # branch 7: node 5 -> node 6
    A[8, 7] = 1                         # branch 8: -> node 7
    A[9, 5], A[9, 7] = 1, -1            # branch 9: node 5 -> node 7
    A[10, 6] = 1                        # branch 10: -> node 6
    A[11, 6] = 1                        # branch 11: -> node 6

    # G: conductance matrix
    G = [G['cv wall']['out'],           # branch 0
         2 * G['cd']['concrete'],       # branch 1
         2 * G['cd']['concrete'],       # branch 2
         2 * G['cd']['insulation'],     # branch 3
         2 * G['cd']['insulation'],     # branch 4
         G['LW'],                       # branch 5
         G['cv wall']['in'],            # branch 6
         G['cv glass']['in'],           # branch 7
         G['cdv glass'],                # branch 8
         2 * G['cd']['glass'],          # branch 9
         G['adv'],                      # branch 10
         G['Kp']]                       # branch 11
    G = np.diag(G)

    # C: capacity matrix
    C = [0,                             # node 0
         C['concrete'],                 # node 1
         0,                             # node 2
         C['insulation'],               # node 3
         0,                             # node 4
         0,                             # node 5
         C['air'],                      # node 6
         C['glass']]                    # node 7
    C = np.diag(C)

    # Sources
    b = np.zeros(nq)                    # all branches
    b[[0, 8, 10, 11]] = 1               # branches with temperature sources

    f = np.zeros(nθ)                    # all nodes
    f[[0, 4, 6, 7]] = 1                 # nodes with flow sources

    # Output temperatures
    y = np.zeros(nθ)                    # all nodes
    y[[6]] = 1                          # temperature nodes of interest

    return A, G, b, C, f, y


def step_response(A, B, C, D, To, Tisp, method='exp'):
    """
    Step response of a state-space model
    for a duration of about 3 - 4 the largest time constant.

    Time step: Δt ≤ 2 * min(-1 / λ), λ - enigenvalues of state matrix A.

    Time duration: t_response ≥ 3 * max(-1 / λ)

    Function calls:
        - time_step()
        - time_response()
        - <integration_Euler()

    Parameters
    ----------
    A : np.Array of float
        State matrix.
    B : np.Array of float
        Input matrix.
    C : np.Array of float
        Output matrix.
    D : np.Array of float
        Feedthrough matrix.
    To : int or float
        outdoor temperature, °C.
    Tisp : int or float, °C
        indoor temperature set point.
    method : str
        'imp': Euler implicit method;
        'exp': Euler explicit method (by default).

    Returns
    -------
    t : Array of floats
        Time vector, seconds.
    y : Array of floats
        Output temperature(s), °C, at time steps t, seconds.
        y.shape[0] -> node; y.shape[1] -> time
    """
    λ = np.linalg.eig(A)[0]             # eigenvalues of matrix A

    # Time step
    dtmax = 2 * min(-1. / λ)            # maximum time step for stability
    dt = time_step(0.5 * dtmax)         # rounded time step

    # Time duration
    settling_time = 3 * max(-1. / λ)
    duration = time_response(settling_time)

    n = int(np.floor(duration / dt))    # number of time steps
    t = np.arange(0, n * dt, dt)        # time vector for n time steps

    # Input vector
    u = np.zeros([np.shape(D)[1], n])   # u = [To To To Tisp Φo Φi Qa Φa]
    u[0:3, :] = To * np.ones([3, n])    # To for n time steps
    u[3, :] = Tisp * np.ones([1, n])    # Tisp for n time steps

    y = integration_Euler(A, B, C, D, t, dt, u, x0=0, method=method)

    yt = pd.Series(y.T[:, 0], index=pd.to_datetime(t, unit='s'))

    return yt


def step_plot(G, C):
    """
    Plots step response.

    Calls :
        - thermal_circuit()
        - dm4bem.tc2ss()
        - <step_response():
            - time_step()
            - time_response()
            - <integration_Euler()

    Parameters
    ----------
    G : dict
        Conductances given by function conductance_capacity()
        Thermal conductances, W/K :
         - ['cd']:
             - ['concrete']
             - ['insulation']
             - ['glass']
         - ['cv wall']:
             - ['in']
             - ['out'];
         - ['cv glass']:
             - ['in']
             - ['out'];
         - ['LW']
         - ['adv']
         - ['cdv glass']

    C : pandas.Series
        Capacities given by function conductance_capacity(), J/kg:
            - ['concrete']
            - ['insulation']
            - ['air']

    Returns
    -------
    None.

    """
    A, G, b, C, f, y = thermal_circuit(G, C)
    As, Bs, Cs, Ds = dm4bem.tc2ss(A, G, b, C, f, y)

    start = time.time()
    y_exp = step_response(As, Bs, Cs, Ds, To=10, Tisp=20, method='exp')
    y_imp = step_response(As, Bs, Cs, Ds, To=10, Tisp=20, method='imp')
    end = time.time()

    df = pd.concat([y_exp, y_imp], axis=1)
    df.columns = ['Explicit', 'Implcit']
    df.plot()

    print(f'Time step = {df.index[1] - df.index[0]} \t'
          f'execution time = {end - start:.3f} s')


def weather_response(A, B, C, D,
                     start_date, end_date,
                     filename,
                     surface_orientation,
                     albedo,
                     Tisp, Ti0, Qa,
                     method):
    """
    Reponse of state-space model to outdoor temperature and solar radiation
    from weather file.

    Calls:
        - time_step()
        - dm4bem.read_epw()
        - dm4bem.sol_rad_tilt_surf()
        - <integration_Euler()

    Parameters
    ----------
    A : Array of float
        State matrix.
    B : Array of float
        Input matrix.
    C : Array of float
        Output matrix.
    D : Array of float
        Feedthrough matrix.
    start_date : datetimes, YYYY-MM-DD HH:MM:SS
        Start date and time of the simulation.
    end_date : datetimes, YYYY-MM-DD HH:MM:SS
        End date and time of the simulation.
    filename : string
        Name of the `.epw` EnergyPlus weather file with folder path.
        Used by dm4bem.read_epw().
    surface_orientation : dict{'slope', 'azimuth', 'latitude'}
        Orientation of the surface for calculating solar radiation
        (diffuse, direct, reflected). Used by dm4bem.sol_rad_tilt_surf().
    albedo : float, adimensional [0, 1]
        Reflection coefficnet of diffuse solar radiation.
    Tisp : float
        Indoor temperature set-point, °C.
        In this implementation, Tisp is constant.
    Ti0 : float
        Initial value of the indoor temperature for time integration, °C.
    Qa : float
        Auxiliary heat flows in the indoor space, W.
        In this implementation, Qa is constant.
    method : str
        Type of Euler integration: 'explicit' or 'implicit'.
        Only the 1st letter ('e' or 'i') counts.

    Returns
    -------
    t : numpy.ndarray
        Time vector, seconds.
    y : numpy.ndarray
        Output of the model, indoor temperature, °C, at time steps t.
    q_HVAC : pandas.Series
        Time series of the heat load, W.
        q_HVAC > 0 -> heating
        q_HVAC < 0 -> cooling
    To : pandas.Series
        Time series of the outdoor temperature, °C, from the weather file.
    Etot : pandas.Series
        Time series of total solar radiation, W/m², from the weather file.
    """
    λ = np.linalg.eig(A)[0]             # eigenvalues of matrix A

    # Time step
    dtmax = 2 * min(-1. / λ)                # maximum time step for stability
    dt = time_step(0.5 * dtmax)         # rouded time step

    # Input data
    # weather data
    [data, meta] = dm4bem.read_epw(filename, coerce_year=None)
    weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
    del data

    weather.index = weather.index.map(lambda t: t.replace(year=2000))
    weather = weather[(weather.index >= start_date) & (
        weather.index < end_date)]

    rad_surf = dm4bem.sol_rad_tilt_surf(
        weather, surface_orientation, albedo)
    rad_surf['Etot'] = rad_surf.sum(axis=1)

    # resample weather data
    data = pd.concat([weather['temp_air'], rad_surf['Etot']], axis=1)
    data = data.resample(str(dt) + 'S'
                         ).bfill(limit=1).interpolate(method='linear')
    data = data.rename(columns={'temp_air': 'To'})

    # indoor temperature set point and auxiliary heat
    data['Ti'] = Tisp * np.ones(data.shape[0])
    data['Qa'] = Qa * np.ones(data.shape[0])

    # input vector
    To = data['To']
    Ti = data['Ti']
    Φo = αS['wall out'] * data['Etot']
    Φi = αS['wall in'] * data['Etot']
    Qa = data['Qa']
    Φa = αS['glass'] * data['Etot']

    Φo = Φo.rename('Φo')

    u = pd.concat([To, To, To, Ti, Φo, Φi, Qa, Φa], axis=1)
    u = u.to_numpy().T

    # Time integration
    t = dt * np.arange(data.shape[0])           # time vector
    y = integration_Euler(A, B, C, D, t, dt, u, x0=Ti0, method=method)

    q_HVAC = G['Kp'] * (data['Ti'] - y[0, :])   # heat load, W
    q_HVAC = q_HVAC.rename('q_HVAC')
    # Etot = data['Etot']                       # total solar iradiance, W/m²

    # return t, y, q_HVAC, To, Φo
    df = pd.concat([q_HVAC, To, Φo], axis=1)
    df['y'] = y.T
    return df


def weather_plot(G, C,
                 start_date, end_date,
                 filename,
                 surface_orientation,
                 albedo=0.2,
                 Tisp=20, Ti0=20, Qa=0,
                 method='explicit'):
    """
    Plots simulation results obtained by `weather_response()`.
    Gives two plots fnction of time:
        1. Indoor temperature and outdoor temperature, °C.
        2. Total solar radiation absorbed and heat load, W.

    Function structure:
        - thermal_circuit()
        - dm4bem.tc2ss()
        - weather_response():
            - time_step()
            - dm4bem.read_epw()
            - dm4bem.sol_rad_tilt_surf()
            - integration_Euler():

    Parameters
    ----------
    G : pandas.Series
        Conductances given by function conductance_capacity()
        Thermal conductances, W/K :
         - ['cd']:
             - ['concrete']
             - ['insulation']
             - ['glass']
         - ['cv wall']:
             - ['in']
             - ['out'];
         - ['cv glass']:
             - ['in']
             - ['out'];
         - ['LW']
         - ['adv']
         - ['cdv glass']
    C : pandas.Series
        Capacities given by function conductance_capacity(), J/kg:
            - ['concrete']
            - ['insulation']
            - ['air']
    start_date : datetimes, YYYY-MM-DD HH:MM:SS
        Start date and time of the simulation.
    end_date : datetimes, YYYY-MM-DD HH:MM:SS
        End date and time of the simulation.
    filename : string
        Name of the `.epw` EnergyPlus weather file with folder path.
        Used by dm4bem.read_epw().
    surface_orientation : dict{'slope', 'azimuth', 'latitude'}
        Orientation of the surface for calculating solar radiation
        (diffuse, direct, reflected). Used by dm4bem.sol_rad_tilt_surf().
    albedo : float, adimensional [0, 1]
        Reflection coefficnet of diffuse solar radiation.
    Tisp : float
        Indoor temperature set-point, °C.
        In this implementation, Tisp is constant.
    Ti0 : float
        Initial value of the indoor temperature for time integration, °C.
    Qa : float
        Auxiliary heat flows in the indoor space, W.
        In this implementation, Qa is constant.
    method : str
        Type of Euler integration: 'explicit' or 'implicit'.
        Only the 1st letter ('e' or 'i') counts.

    Returns
    -------
    None.

    """

    A, G, b, C, f, y = thermal_circuit(G, C)
    As, Bs, Cs, Ds = dm4bem.tc2ss(A, G, b, C, f, y)

    df = weather_response(As, Bs, Cs, Ds,
                          start_date, end_date,
                          filename,
                          surface_orientation,
                          albedo,
                          Tisp, Ti0, Qa,
                          method)
    pass
    fig, axes = plt.subplots(nrows=2, ncols=1)

    df[['y', 'To']].plot(ax=axes[0],
                         ylabel='Temperarure, °C')

    df[['q_HVAC', 'Φo']].plot(ax=axes[1],
                              xlabel='Time',
                              ylabel='Power, W')
    plt.show()


def integration_Euler(A, B, C, D, t, dt, u, x0=0, method='exp'):
    """

    Parameters
    ----------
    A : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    C : TYPE
        DESCRIPTION.
    D : TYPE
        DESCRIPTION.
    duration : TYPE
        DESCRIPTION.
    dt : TYPE
        DESCRIPTION.
    u : TYPE
        DESCRIPTION.
    x0 : TYPE
        DESCRIPTION.
    method : TYPE
        DESCRIPTION.

    Returns
    -------
    y : TYPE
        DESCRIPTION.

    """

    # Numerical integration in time

    n_s = A.shape[0]                        # number of state variables
    x = np.zeros([n_s, t.shape[0]])         # x in time
    x[:, 0] = x0                            # initial condition

    I = np.eye(n_s)                         # identity matrix

    for k in range(np.shape(t)[0] - 1):
        if method[0] == 'expplicit'[0]:        # Euler explicit
            x[:, k + 1] = (I + dt * A) @\
                x[:, k] + dt * B @ u[:, k]
        else:                               # Euler implicit
            x[:, k + 1] = np.linalg.inv(I - dt * A) @\
                (x[:, k] + dt * B @ u[:, k])
    # Outputs
    y = C @ x + D @  u
    return y


def time_step(dtmax):
    """

    Parameters
    ----------
    dtmax : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if dtmax < 1:
        print('Exit: dt_max < 1 s')
        sys.exit()
    elif dtmax < 10:                        # round to 1 sec
        dt = np.floor(dtmax)
    elif dtmax < 60:                        # round at 10 sec
        dt = 10 * np.floor(dtmax / 10)
    elif dtmax < 360:
        dt = 60 * np.floor(dtmax / 60)      # round at min
    elif dtmax < 3600:
        dt = 360 * np.floor(dtmax / 360)    # round at 10 min
    else:
        dt = 3600 * np.floor(dtmax / 3600)  # round at 1 h

    return dt


def time_response(tmin):
    """

    Parameters
    ----------
    dtmax : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if tmin < 1:
        print('dt_max < 1 s')
        exit()
    elif tmin < 60:                        # round at sec
        dt = np.ceil(tmin)
    elif tmin < 360:
        dt = 60 * np.ceil(tmin / 60)      # round at min
    elif tmin < 3600:
        dt = 360 * np.ceil(tmin / 360)    # round at 10 min
    else:
        dt = 3600 * np.ceil(tmin / 3600)  # round at 1 h

    return dt


# Main
start_time_all = time.time()

G, C, mesh, αS = conductance_capacity(Kp=0)
step_plot(G, C)

# G, C, mesh, αS = conductance_capacity(Kp=1e3)
# step_plot(G, C)

G, C, mesh, αS = conductance_capacity(Kp=0)
start_time = time.time()
weather_plot(G, C,
             start_date='2000-01-03 12:00:00',
             end_date='2000-02-05 18:00:00',
             filename='../weather_data/FRA_Lyon.074810_IWEC.epw',
             surface_orientation={'slope': 90,
                                  'azimuth': 0,
                                  'latitude': 45},
             albedo=0.2,
             Tisp=20, Ti0=20, Qa=0,
             method='explicit')
print(f"Weather, Kp = 0: time = {time.time() - start_time:.2f} s")

G, C, mesh, αS = conductance_capacity(Kp=1e2)
start_time = time.time()
weather_plot(G, C,
             start_date='2000-01-03 12:00:00',
             end_date='2000-02-05 18:00:00',
             filename='../weather_data/FRA_Lyon.074810_IWEC.epw',
             surface_orientation={'slope': 90,
                                  'azimuth': 0,
                                  'latitude': 45},
             albedo=0.2,
             Tisp=20, Ti0=20, Qa=0,
             method='implicit')
print(f"Weather, Kp = 1e2: time = {time.time() - start_time:.2f} s")

# G, C, mesh, αS = conductance_capacity(Kp=1e4)
# start_time = time.time()
# weather_plot(G, C,
#              start_date='2000-01-03 12:00:00',
#              end_date='2000-02-05 18:00:00',
#              filename='../weather_data/FRA_Lyon.074810_IWEC.epw',
#              surface_orientation={'slope': 90,
#                                   'azimuth': 0,
#                                   'latitude': 45},
#              albedo=0.2,
#              Tisp=20, Ti0=20, Qa=0,
#              method='implicit')
# print(f"Weather w/ air c, Kp = 1e4: time = {time.time() - start_time:.2f} s")

G, C, mesh, αS = conductance_capacity(Kp=1e4,
                                      air={'ρ': 1.2,
                                           'c': 0},)
start_time = time.time()
weather_plot(G, C,
             start_date='2000-01-03 12:00:00',
             end_date='2000-02-05 18:00:00',
             filename='../weather_data/FRA_Lyon.074810_IWEC.epw',
             surface_orientation={'slope': 90,
                                  'azimuth': 0,
                                  'latitude': 45},
             albedo=0.2,
             Tisp=20, Ti0=20, Qa=0,
             method='implicit')
print(f"Weather w/o air c, Kp = 1e4: time = {time.time() - start_time:.2f} s")

print(f"Time all = {time.time() - start_time_all:.2f} s")