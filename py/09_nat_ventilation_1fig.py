#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 08:45:51 2023

@author: cghiaus

Estimate the flow rate by ventilation

References
ASHARE Fundamentals (2009) Chapter 16 Ventilation and infiltration SI_F9_Ch16
Solve one nonliear equation
https://faculty.math.illinois.edu/~hirani/cbmg/nonlinear.html
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# Data
# standard acceleration due to gravity, m/s²
g = 9.8
# pressure exponent, dimensionless
no, nc = 0.64, 0.50             # orrifice, cheminée
# flow coefficients, m³/(s·Paⁿ)
Ko, Kc = 10 / 3600, 50 / 3600   # orrifice, chimney
# density, kg/m³
ρi, ρo = 1.191, 1.315           # indoor, outdoor
# pressure coefficients, dimensioneless
Cpn, Cps, Cpc = 0.4, -0.55, -0.7    # north, south, chimney
# hignt, m
z = 8

# Static pressure, Pa
Psi = -ρi * g * z   # indoors
Pso = -ρo * g * z   # outdoors


def f(pi):
    """
    Error in mass balance equation as a function of indoor pressure.
    Used to find indoor pressure pi which makes the mass balance zero.
    The north and south flows are enteriing.

    Parameters
    ----------
    pi : float
        Indoor pressure, Pa.

    Returns
    -------
    y : float
        Mass flow rate unbalanced in the mass balance equation
    """
    y = 2 * ρo * Ko * (Pdn - pi)**no \
        + 2 * ρo * Ko * (Pds - pi)**no \
        - ρi * Kc * (pi + Psi - Pdc - Pso)**nc
    return y


# Without wind
# ================================================
print('Without wind')
v = 0       # m/s, wind speed

# dynamic pressures, Pa
Pdn = Cpn * 0.5 * ρo * v**2     # north face
Pds = Cps * 0.5 * ρo * v**2     # south face
Pdc = Cpc * 0.5 * ρo * v**2     # chimney

pi_min, pi_max = -(Psi - Pdc - Pso), min(Pdn, Pds)
print(f"Pressure domain: {pi_min:.2f} Pa < pi < {pi_max:.2f} Pa")

# plot error in mass balance
pi = np.linspace(pi_min, pi_max, 100)

fig, ax = plt.subplots()

without_wind, = ax.plot(pi, f(pi), 'b')
ax.axvline(x=pi_min, color='b', **{'linestyle': 'dashed'})
ax.axvline(x=pi_max, color='b', **{'linestyle': 'dashed'})

# numerical solution
root = float(fsolve(f, np.mean([pi_min, pi_max])))
print(f"pi = {root:.2f} Pa")

# verification of mass balance, kg/s
mn = ρo * Ko * (Pdn - root)**no
ms = ρo * Ko * (Pds - root)**no
mc = ρi * Kc * (root + Psi - Pdc - Pso)**nc
print("Verify mass balance:")
print(f"mass flow: north: {2 * mn:.2f} kg/s, south: {2 * ms:.2f} kg/s")
print(f"mass flow chimney: {mc:.2f} kg/s")
print(f"error in mass flow balance: {2 * mn + 2 * ms - mc:.2f} kg/s\n")

# With wind
# ================================================
print('With wind')
v = 18 * 1000 / 3600            # m/s, wind speed

# dynamic pressures, Pa
Pdn = Cpn * 0.5 * ρo * v**2     # north face
Pds = Cps * 0.5 * ρo * v**2     # south face
Pdc = Cpc * 0.5 * ρo * v**2     # chimney

pi_min, pi_max = -(Psi - Pdc - Pso), min(Pdn, Pds)
print(f"Pressure domain with wind: {pi_min:.2f} Pa < pi < {pi_max:.2f} Pa")

# plot error in mass balance
pi = np.linspace(pi_min, pi_max, 100)

with_wind, = ax.plot(pi, f(pi), 'g')
ax.axhline(color='k')
ax.axvline(x=pi_min, color='g', **{'linestyle': 'dashed'})
ax.axvline(x=pi_max, color='g', **{'linestyle': 'dashed'})
ax.set_ylabel(r'Error in mass balance $\dot{m}$ [kg/s]')
ax.set_xlabel(r'Indoor pressure $p_i$ [Pa]')

ax.axhline(color='k')
ax.set_ylabel(r'Error in mass balance $\dot{m}$ [kg/s]')
ax.set_xlabel(r'Indoor pressure $p_i$ [Pa]')
ax.grid()
ax.legend((without_wind, with_wind), ('Without wind', 'With wind'))
# fig.tight_layout()

# numerical solution
root = float(fsolve(f, np.mean([pi_min, pi_max])))
print(f"pi = {root:.2f} Pa")

mn = ρo * Ko * (Pdn - root)**no
ms = ρo * Ko * (Pds - root)**no
mc = ρi * Kc * (root + Psi - Pdc - Pso)**nc
print("Verify mass balance:")
print(f"mass flow: north: {2 * mn:.2f} kg/s, south: {2 * ms:.2f} kg/s")
print(f"mass flow chimney: {mc:.2f} kg/s")
print(f"error in mass flow balance: {2 * mn + 2 * ms - mc:.2f} kg/s\n")
