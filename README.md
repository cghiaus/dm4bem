# Dynamic Models for Building Energy Management

These tutorials present a complete example of thermal dynamic simulation of a building. For the sake of simplicity, a [toy model](./figures/03_cube_principle.png) is used for the building in which 5 identical two-layer walls, a glass wall, air infiltration, and an indoor temperature control system are modelled.

The tutorials go through obtaining weather data from internet, modelling the thermal transfer with thermal networks, transforming the thermal networks into systems of differential algebraic equations and state-space representation, and implementing customized control algorithms into the numerical integration loop.

The advantage of the method, as compared with other existing alternatives, is that the state-space representation is obtained; therefore eigenvalues analysis is achievable.

The disadvantage is that, in the current implementation, application on large models is tedious and prone to errors. 

The notebooks can be run interactively on `MyBinder.org` by clicking on the button below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/dm4bem/HEAD)

## Tutorials
- 01: [Weather data and solar radiation on a tilted surface](01WeatherData.ipynb).
- 02: [Thermal circuit and state-space representation for a thermal circuit with capacities in every node: simple wall](02SimpleWall.ipynb).
- 03: [Thermal circuit and state-space representation for a thermal circuit with capacities in some nodes: cubic building](03CubicBuilding.ipynb).
- 04: [Thermal circuits assembling](04AssemblingCircuits.ipynb).
- 05: [Switch between models: heating & cooling and free-running](05SwitchModels.ipynb).
- 06: [Control input: heating & cooling and free-running](06Control_Input.ipynb).

## Sessions
1. **Model**
 - Draw the plan of a simple building.
 - Formulate the hypothesis for boundary conditions, windows, doors, and wall composition.
 - Write down the incidence matrix $A$, the conductance matrix $G$ and the capacity matrix $C$.
 - Define the inputs: temperature sources (vector **b**) and flow rate sources (vector **f**).
2. **Pyhton implementation: steady-state**
 - Implement the matrices **A**, **G** and **C**.
 - Calculate the solar flows.
 - Write the input vectors **b** and **f** in time.
 - Calculate steady-state response.
3. **Pyhton implementation: simulation**
 - Simulate a step response.
 - Simulate the response to weather data.
 - Debug and optimize.
 - Complex controllers (dead-band, model predictive control).
4. **Write reproducible report**
 - Write the report in *Jupyter* notebooks.
 - Publish the report on *GitHub* and *MyBinder*.


**References**

1. C. Ghiaus (2013). Causality issue in the heat balance method for calculating the design heating and cooling load. *Energy* 50: 292-301
[DOI 10.1016/j.energy.2012.10.024](http://dx.doi.org/10.1016/j.energy.2012.10.024), [HAL 03605823]( https://hal.archives-ouvertes.fr/hal-03605823/document)

2. C. Ghiaus, N. Ahmad (2020). Thermal circuits assembling and state-space extraction for modelling heat transfer in buildings, *Energy*, 195:117019
[DOI 10.1016/j.energy.2020.117019](https://doi.org/10.1016/j.energy.2020.117019), [HAL 03600778](https://hal.archives-ouvertes.fr/hal-03600778/document)

3. C. Ghiaus (2021). Dynamic Models for Energy Control of Smart Homes, in *S. Ploix M. Amayri, N. Bouguila (eds.) Towards Energy Smart Homes*, Online ISBN: 978-3-030-76477-7, Print ISBN: 978-3-030-76476-0, Springer, pp. 163-198 (ref.)
[DOI 10.1007/978-3-030-76477-7_5](https://doi.org/10.1007/978-3-030-76477-7_5), [HAL 03578578](https://hal.archives-ouvertes.fr/hal-03578578/document)

4. J. Kneifel (2013). Annual Whole Building Energy Simulation of the NIST Net Zero Energy Residential Test Facility Design, *NIST Technical Note 1767*, [DOI 10.6028/NIST.TN.1767](https://doi.org/10.6028/NIST.TN.1767)

5. U.S. Department of Energy (2022). EnergyPlus v.22.1.0 Documentation, Engineering Reference ([link](https://energyplus.net/assets/nrel_custom/pdfs/pdfs_v22.1.0/EngineeringReference.pdf))

# Exam questions
1. Defintion of science.
2. Shortly describe the four paradigms of science (empirical, theoretical, computational and data science).
3. Reproducibility crises: definition and how to overcome it.
4. Describe the reproducibility spectrum in Computational Science.
5. Definition of physical and computational causality.
6. Definition of inputs, outputs and states of a dynamic model.
7. Conservation laws: two examples.
8. Relation between conservation laws and symmetry in physics.
9. Constitutive laws: one example of a universal law and one of a phenomenological law.
10. Explain why there are only seven fundamental units in the SI system of units.
11. What is the difference between the classical SI system of units and the system addopted on 20 May 2019? Why is this difference important?
12. What is the relationship between energy and temperature?
13. Draw the basic network for heat transfer modelling. Explain each element of the network:
    - temperaure nodes,
    - flow branches,
    - conductances,
    - capacities,
    - temperature sources,
    - flow sources.
14. Draw the framework for obtaining the *difussion equation*.
15. Show the analogy between:
    - heat transfer,
    - mass transfer,
    - electrical conduction.
16. The framework for obtaining models for transport phenomena (i.e.n heat trasfer, mass transfer, momentum transfer, electrical conduction).
17. Define the modes of heat transfer and give the expression of conductance for:
    - conduction,
    - convection,
    - radiation,
    - advection.
18. Conservation of energy in steady-state and in dynamics.
19. Definition of sensible heat.
29. Surface phenomena and volume phenomena in energy balance equation.
21. Draw a wall and a window. Make a thermal network model of this system.
22. Explain the difference between *Differential Algebraic Equations* model and *state-space* representation.
23. Explain the shape and the elements of the following matrices and vectors of a *Differential Algebraic Equations* model $C \dot{\theta} = -A^T G A \theta + A^T G b + f$:
    - $A$ - incidence marix;
    - $G$ - conductance matrix;
    - $C$- capacity matrix;
    - $b$ - temperature source vector;
    - $f$ - flow-rate source vector.

# Written report
The report will be written in *Jupyter* notebook, posted on *GitHub.com* and liked to *MyBinder.org*.

The general structure of the report:
- Front page: title, author(s), date.
- Contents.
- Description of the building: drawing, dimensions, materials, material properties, etc.
- Hypothesis: location, boundary conditions, schedule for usage, etc.
- Thermal model (with justifications).
- Mathematical model: Differential Algebraic Equations (matrices $A$, $G$, and $C$, vectors $b$ and $f$) and state-space representation (matrices $A$, $B$, $C$ and $D$ and vector $u$).
- Model implementation in Python.
- Steady-state results.
- Dynamic simulation results.
- Optimization (e.g., insulation, HVAC control, ventilation rate, solar shading).


**Licence**

Code is released under [MIT Lincence](https://choosealicense.com/licenses/mit/).

[![Creative Commons License](http://i.creativecommons.org/l/by/4.0/88x31.png)](http://creativecommons.org/licenses/by/4.0/)