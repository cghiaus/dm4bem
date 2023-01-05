# Tutorials for MATLAB®/Octave scripts

The scripts are written for MATLAB® v.6 or GNU Octave v. 6.

1. `t01ReadWeather.pdf` *Weather data and solar radiation*: obtain weather data from Energy Plus weather files and find solar radiation on a tilted surface.
2. `t02SimpleWall.pdf`  *Simple wall*: Differential-algebraic equations (DAE) model of a two layer wall with Dirichelt end von Neumann boundary conditions.
3. `t03Cube2wFB.pdf` *Cube with two different walls: feed-back indoor temperature control*: model of a cubic building with 5 identical walls and one glass wall. HVAC controller, air infiltration and internal sources are modeled.
4. `t04Cube2wFBmesh.pdf` *Cube with two different walls with variable meshgrid*: similar to `t03Cube2wFB.pdf`with the difference that the meshing of the wall layers is variable.
5. `t05CubeFBAss.pd` *Assembling thermal circuits*: similar to `t03Cube2wFB.pdf` with the difference that the thermal model is obtained by assembling thermal networks.
6. `t06CubeFBHeat.pdf` *Cubic room heated by a fan-coil*: similar to `t03Cube2wFB.pdf`with the difference that the controler has a dead-band (non-linear controller).