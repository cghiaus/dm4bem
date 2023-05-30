# Dynamic Models for Building Energy Management

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/dm4bem/HEAD)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/cghiaus/dm4bem/blob/main/LICENSE)


These tutorials present a complete example of thermal dynamic simulation of a building. For the sake of simplicity, a [toy model](https://en.m.wikipedia.org/wiki/Toy_model) is used for the [building](./figures/03_cube_principle.png) in which 5 identical two-layer walls, a glass wall, air infiltration, and an indoor temperature control system are modelled.

The tutorials go through obtaining weather data from internet, modelling the thermal transfer with thermal networks, transforming the thermal networks into systems of differential algebraic equations and state-space representation, and implementing customized control algorithms into the numerical integration loop.

The advantage of the method, as compared with other existing alternatives, is that the state-space representation is obtained; therefore eigenvalues analysis is achievable.

The disadvantage is that, in the current implementation, application on large models is tedious and prone to errors. 

The notebooks can be run interactively on `MyBinder.org` by clicking on the button below:

*Note*: The repository containes the MATLAB® version 6 scripts and the PDFs of the tutorials in the folder `M`.

## Tutorials
1. [Weather data and solar radiation on a tilted surface](01WeatherData.ipynb).
2. [Thermal circuit and state-space representation for a thermal circuit with capacities in every node: simple wall](02SimpleWall.ipynb).
3. [Thermal circuit and state-space representation for a thermal circuit with capacities in some nodes: cubic building](03CubicBuilding.ipynb).
4. [Thermal circuits assembling](04AssemblingCircuits.ipynb).
5. [Switch between models: heating & cooling and free-running](05SwitchModels.ipynb).
6. [Control input: heating & cooling and free-running](06Control_Input.ipynb).
7. [Radiation coupled with convection](07Coupled_rad_convection.ipynb).
8. [Sensible thermal load in steady-state](08Thermal_load.ipynb)
9. [Air flow by ventilation](09Air_flow_ventilation.ipynb)

[Slides of the lectures](https://filesender.renater.fr/?s=download&token=e4432033-c204-4d44-b23a-e47d3265012d)

## Project sessions
The assignments are written in *Jupyter* notebooks and posted on Github in the repository indicated in each assignment. Each repository needs to have:
- `README.md` file that contains at least the names of the members of the team and a link to `mybinder.org`.
- `environment.yml` file that lists the dependencies required for the project.

The *Jupyter* notebooks need to contain [Markdown cells](https://jupyter-notebook.readthedocs.io/en/stable/examples/Notebook/Working%20With%20Markdown%20Cells.html) structured with [headings](https://www.markdownguide.org/basic-syntax/#headings) and equations need to be written in [LaTeX](https://jupyter-notebook.readthedocs.io/en/stable/examples/Notebook/Working%20With%20Markdown%20Cells.html#LaTeX-equations).


1. **Model**
 - Draw the plan of a two-zone building.
 - Formulate the hypothesis for boundary conditions.
 - Chose the types of windows, doors, and walls.
 - Draw the thermal circuit:
     - temperature nodes,
     - flow-rate paths,
     - thermal conductances for conduction, convection, long-wave radiation, advection, and P-controllers,
     - sources of temperature and flow-rate,
- Number the temperature nodes and the flow-rate branches (starting from 0).
 - Calculate the thermal conductances for conduction, convection, long-wave radiation, and advection.
 - Calculate the thermal capacities.
 - Write down the incidence matrix $A$, the conductance matrix $G$ and the capacity matrix $C$ of the system of Algebraic Differential Equations (DAE).
 - Define the inputs: temperature sources (vector $b$) and flow rate sources (vector $f$).
 - Write in Pyhthon the incidence matrix $A$, the conductance matrix $G$ and the capacity matrix $C$ of the system of Algebraic Differential Equations (DAE).
 - Write in Pyhthon the vectors of pointers to the temperature sources $b$, flow-rate sources $f$, and outputs $y$.
 - **Assignment 1**: Model
     - [Grup 1, 2](https://classroom.github.com/a/Fh4jnCT2)
     - [Grup 3, 4](https://classroom.github.com/a/bl3pb-J7)
 
2. **Steady-state**
 - Implement in Python the matrices $A$, $G$ and $C$ of the system of Diferential Algebraic Equations (DAE).
 - Give values to inputs (temperature sources, $b$, and flow rate sources $f$).
 - Calculate steady-state response of the system of Diferential Algebraic Equations (DAE).
 - From the systems of Diferential Algebraic Equations (DAE), obtain the matrices $A_s$, $B_s$, $C_s$, and $D_s$ of the state-space representation.
 - Give the values of the input vector $u = [b_T^T, f_Q^T]^T$.
 - Obtain the steady-state response of the state-space representation.
 - Compare the results obtained for the system of Diferential Algebraic Equations (DAE) with the results obtained for the state-space representation. 
 - **Assignment 2**: Steady-state
     - [Grup 1, 2](https://classroom.github.com/a/6HWx5wze)
     - [Grup 3, 4](https://classroom.github.com/a/T6cSvhT4)
 
3. **Simulate step response**
 - Determine the time step and the settling time.
 - Give the input vector $u$.
 - Integrate in time the differential equations.
 - Plot the results.
 - Discuss the results.
 - **Assignment 3**: Simulate step response
     - [Grup 1, 2](https://classroom.github.com/a/1YlSy6uy)
     - [Grup 3, 4](https://classroom.github.com/a/XR32Fbwz)

4. **Simulate response to weather**
 - Define start and end time.
 - Prepare the inputs:
     - read weather data,
     - calculate the solar irradiance on the walls,
     - resample the weather data
     - give the other inputs (e.g., internal gains),
     - give the input vector in time.
 - Define the initial conditions.
 - Integrate in time.
 - Plot the results.
 - Discuss the results.
 - Implement other controllers (dead-band, model predictive control).
 - Discuss the results.
 - **Assignment 4**: Simulate response to weather
     - [Grup 1, 2](https://classroom.github.com/a/1VJkM4fh)
     - [Grup 3, 4](https://classroom.github.com/a/1VJkM4fh)
 
5. **Reproducible report**
 - Write the report in *Jupyter* notebooks.
 - Publish the report on *GitHub* and *MyBinder*.
 - **Assignment 5**: Reproducible report
     - [Grup 1, 2](https://classroom.github.com/a/4YDPKTYq)
     - [Grup 3, 4](https://classroom.github.com/a/Be3bPaux)

**Support**

1. [GitHub Docs](https://docs.github.com/en):
    - [Getting started with your GitHub account](https://docs.github.com/en/get-started/onboarding/getting-started-with-your-github-account)
    - [Accepting the assignement](https://www.youtube.com/watch?v=jXpT8eOzzCM)
    - [How to add files to Github Repository](https://www.youtube.com/watch?v=-rcz_CZe2Ec)
2. [Anaconda cheetsheet](https://docs.continuum.io/anaconda/user-guide/cheatsheet/)
3. Python
    - [Python tutorial](https://docs.python.org/3/tutorial/#the-python-tutorial)
    - [Python course CS50 (Harvard)](https://cs50.harvard.edu/python/2022/weeks/)
4. [Jupyter notebook cheatsheet](https://medium.com/ibm-data-science-experience/markdown-for-jupyter-notebooks-cheatsheet-386c05aeebed)
5. [NumPy for MATLAB users](http://mathesaurus.sourceforge.net/matlab-numpy.html)
6. [Get started with Binder](https://mybinder.readthedocs.io/en/latest/introduction.html)
7. [Markdown and LaTeX introduction](https://ashki23.github.io/markdown-latex.html)
    - [LaTeX equation editor](https://latex.codecogs.com/eqneditor/editor.php)
    - [Markdown table generator](https://www.tablesgenerator.com/markdown_tables#)


## References

1. C. Ghiaus (2013). Causality issue in the heat balance method for calculating the design heating and cooling load. *Energy* 50: 292-301
[DOI 10.1016/j.energy.2012.10.024](http://dx.doi.org/10.1016/j.energy.2012.10.024), [HAL 03605823]( https://hal.archives-ouvertes.fr/hal-03605823/document)

2. C. Ghiaus, N. Ahmad (2020). Thermal circuits assembling and state-space extraction for modelling heat transfer in buildings, *Energy*, 195:117019
[DOI 10.1016/j.energy.2020.117019](https://doi.org/10.1016/j.energy.2020.117019), [HAL 03600778](https://hal.archives-ouvertes.fr/hal-03600778/document)

3. C. Ghiaus (2021). Dynamic Models for Energy Control of Smart Homes, in *S. Ploix M. Amayri, N. Bouguila (eds.) Towards Energy Smart Homes*, Online ISBN: 978-3-030-76477-7, Print ISBN: 978-3-030-76476-0, Springer, pp. 163-198 (ref.)
[DOI 10.1007/978-3-030-76477-7_5](https://doi.org/10.1007/978-3-030-76477-7_5), [HAL 03578578](https://hal.archives-ouvertes.fr/hal-03578578/document)

4. J. Kneifel (2013). Annual Whole Building Energy Simulation of the NIST Net Zero Energy Residential Test Facility Design, *NIST Technical Note 1767*, [DOI 10.6028/NIST.TN.1767](https://doi.org/10.6028/NIST.TN.1767)

5. U.S. Department of Energy (2022). EnergyPlus v.22.1.0 Documentation, Engineering Reference ([link](https://energyplus.net/assets/nrel_custom/pdfs/pdfs_v22.1.0/EngineeringReference.pdf))

6. Solar Energy Laboratory, University of Wisconsin-Madison (2009). TRNSYS 17 Volume 4 Mathematical Reference ([link](https://web.mit.edu/parmstr/Public/TRNSYS/04-MathematicalReference.pdf))

7. NIST (2008) Guide for the Use of the International System of Units (SI) ([link](https://physics.nist.gov/cuu/pdf/sp811.pdf))

8. BIPM (2008) Evaluation of measurement data — Guide to the expression of uncertainty in measurement ([link](https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf/cb0ef43f-baa5-11cf-3f85-4dcd86f77bd6))

9. BIPM (2019) The International System of Units (SI), 9th edition
([link](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9.pdf))


# Exam questions
1. Definition of science.
2. Shortly describe the four paradigms of science (empirical, theoretical, computational and data science).
3. Reproducibility crises: definition and how to overcome it.
4. Describe the reproducibility spectrum in Computational Science.
5. What is a (mathematical) model?
6. Definition of physical and computational causality.
7. Definition of inputs, outputs and states of a dynamic model.
8. Conservation laws: two examples.
9. Relation between conservation laws and symmetry in physics.
10. Constitutive laws: one example of a universal law and one of a phenomenological law.
11. Explain why there are only seven fundamental units in the SI system of units.
12. What is the difference between the classical SI system of units and the system adopted on 20 May 2019? Why is this difference important?
13. What is the relationship between energy and temperature?
14. Draw the basic network for heat transfer modelling. Explain each element of the network:
    - temperature nodes,
    - flow branches,
    - conductances,
    - capacities,
    - temperature sources,
    - flow sources.
15. Show the analogy between:
    - heat transfer,
    - mass transfer,
    - electrical conduction.
16. Draw the framework for obtaining models for transport phenomena (i.e., heat trasfer, mass transfer, momentum transfer, electrical conduction).
17. Define the modes of heat transfer and give the expression of conductance for:
    - conduction,
    - convection,
    - radiation,
    - advection.
18. Conservation of energy in steady-state and in dynamics.
29. Definition of sensible heat.
20. Surface phenomena and volume phenomena in energy balance equation.
21. Draw a wall and a window. Make a thermal network model of this system.
22. Explain the difference between *Differential Algebraic Equations* model and *state-space* representation.
23. Explain the shape and the elements of the following matrices and vectors of a *Differential Algebraic Equations* model $C \dot{\theta} = -A^T G A \theta + A^T G b + f$:
    - $A$ - incidence marix;
    - $G$ - conductance matrix;
    - $C$ - capacity matrix;
    - $b$ - temperature source vector;
    - $f$ - flow-rate source vector.
24. What is the relationship between the eigenvalues of the state matrix and the time constants of the thermal model?
25. Define the numerical stability. What is the condition for the numerical stability of Euler explicit method of integration in time?
26. Explain the differences between Euler implicit and explicit methods.
27. Comment on the advantages and the disadvantages of Euler implicit and explicit methods.
28. Consider the the model of heat transfert through a simple wall (tutorial 29). Explain the qualitative differences between the step response to outdoor temperature input and the step response to indoor heat flow rate.
29. Discuss the influence of the initial conditions on the dynamic response.
30. How can be estimated the response (or settling) time of a dynamic system? What is the usefulness of the response time in dynamic simulation?
31. How can be estimated the time step for the dynamic simulation? What can be done if the time step is too small as compared to the time step needed for the problem?

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
