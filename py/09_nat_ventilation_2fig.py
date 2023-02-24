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
# flow coefficients, m³/(s·Paⁿ)
Ko, Kc = 10 / 3600, 50 / 3600   # orrifice, chimney
# pressure exponent, dimensionless
no, nc = 0.64, 0.50             # orrifice, cheminée
# density, kg/m³
ρi, ρe = 1.191, 1.315           # indoor, outdoor
# pressure coefficients, dimensioneless
Cpn, Cps, Cpc = 0.4, -0.55, -0.7    # north, south, chimney
# hignt, m
z = 8

# Static pressure, Pa
Psi = -ρi * g * z   # indoor
Pso = -ρe * g * z   # outdoor


def f(Pi):
    """
    Error in mass balance equation as a function of indoor pressure.
    Used to find indoor pressure Pi which makes the mass balance zero.
    The north and south flows are enteriing.

    Parameters
    ----------
    Pi : float
        Indoor pressure, Pa.

    Returns
    -------
    y : float
        Mass flow rate unbalance in the mass balance equation
    """
    y = 2 * ρe * Ko * (Pdn - Pi)**no \
        + 2 * ρe * Ko * (Pds - Pi)**no \
        - ρi * Kc * (Pi + Psi - Pdc - Pso)**nc
    return y


# Without wind
# ================================================
print('Without wind')
v = 0       # m/s, wind speed

# dynamic pressures, Pa
Pdn = Cpn * 0.5 * ρe * v**2     # north face
Pds = Cps * 0.5 * ρe * v**2     # south face
Pdc = Cpc * 0.5 * ρe * v**2     # chimney

Pi_min, Pi_max = -(Psi - Pdc - Pso), min(Pdn, Pds)
print(f"Pressure domain: {Pi_min:.2f} Pa < Pi < {Pi_max:.2f} Pa")

# plot error in mass balance
Pi = np.linspace(Pi_min, Pi_max, 100)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
ax1.axhline(color='k')
ax1.axvline(x=Pi_min, color='r', **{'linestyle': 'dashed'})
ax1.axvline(x=Pi_max, color='r', **{'linestyle': 'dashed'})
ax1.set_ylabel(r'Error in mass balance $\dot{m}$ [kg/s]')
ax1.set_xlabel(r'Indoor pressure $p_i$ [Pa]')
ax1.grid()
ax1.plot(Pi, f(Pi), 'blue')
ax1.set_title('Without wind')

# numerical solution
root = float(fsolve(f, np.mean([Pi_min, Pi_max])))
print(f"Pi = {root:.2f} Pa")

# verification of mass balance, kg/s
mn = ρe * Ko * (Pdn - root)**no
ms = ρe * Ko * (Pds - root)**no
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
Pdn = Cpn * 0.5 * ρe * v**2     # north face
Pds = Cps * 0.5 * ρe * v**2     # south face
Pdc = Cpc * 0.5 * ρe * v**2     # chimney

Pi_min, Pi_max = -(Psi - Pdc - Pso), min(Pdn, Pds)
print(f"Pressure domain with wind: {Pi_min:.2f} Pa < Pi < {Pi_max:.2f} Pa")

# plot error in mass balance
Pi = np.linspace(Pi_min, Pi_max, 100)

ax2.axhline(color='k')
ax2.axvline(x=Pi_min, color='r', **{'linestyle': 'dashed'})
ax2.axvline(x=Pi_max, color='r', **{'linestyle': 'dashed'})
ax2.set_ylabel(r'Error in mass balance $\dot{m}$ [kg/s]')
ax2.set_xlabel(r'Indoor pressure $p_i$ [Pa]')
ax2.grid()
ax2.plot(Pi, f(Pi), 'blue')
ax2.set_title('With wind')

# fig, ax = plt.subplots(nrows=1, ncols=1)
# ax.plot(Pi, f(Pi), 'blue')
# ax.legned('With wind')

# ax.axhline(color='k')
# ax.axvline(x=Pi_min, color='r', **{'linestyle': 'dashed'})
# ax.axvline(x=Pi_max, color='r', **{'linestyle': 'dashed'})
# ax.set_ylabel(r'Error in mass balance $\dot{m}$ [kg/s]')
# ax.set_xlabel(r'Indoor pressure $p_i$ [Pa]')
# ax.grid()

# numerical solution
root = float(fsolve(f, np.mean([Pi_min, Pi_max])))
print(f"Pi = {root:.2f} Pa")

mn = ρe * Ko * (Pdn - root)**no
ms = ρe * Ko * (Pds - root)**no
mc = ρi * Kc * (root + Psi - Pdc - Pso)**nc
print("Verify mass balance:")
print(f"mass flow: north: {2 * mn:.2f} kg/s, south: {2 * ms:.2f} kg/s")
print(f"mass flow chimney: {mc:.2f} kg/s")
print(f"error in mass flow balance: {2 * mn + 2 * ms - mc:.2f} kg/s\n")

fig.tight_layout()
plt.savefig("Figure_1.svg")
