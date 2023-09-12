#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 10:51:12 2022
Last updated on Tue Jul  4 14:42:21 2023

@author: cghiaus

Import functions for EPW data files.
Adapted from
https://github.com/pvlib/pvlib-python/blob/master/pvlib/iotools/epw.py
"""

import numpy as np
import pandas as pd
import os
import glob
import ast

# Physical constants
σ = 5.670_374_419e-8    # W⋅m⁻²⋅K⁻⁴ Stefan-Bolzmann constant


def read_epw(filename, coerce_year=None):
    '''
    Read an Energy Plus Weather (EPW) file into a pandas dataframe.

    Note that values contained in the metadata dictionary are unchanged
    from the EPW file.

    EPW files are commonly used by building simulation professionals
    and are widely available on the web. For example via:
    https://energyplus.net/weather , http://climate.onebuilding.org or
    http://www.ladybug.tools/epwmap/


    Parameters
    ----------
    filename : String
        Can be a relative file path or absolute file path.

    coerce_year : None or int, default None
        If supplied, the year of the data will be set to this value. This can
        be a useful feature because EPW data is composed of data from
        different years.
        Warning: EPW files always have 365*24 = 8760 data rows;
        be careful with the use of leap years.


    Returns
    -------
    data : DataFrame
        A pandas dataframe with the columns described in the table
        below. For more detailed descriptions of each component, please
        consult the EnergyPlus Auxiliary Programs documentation [1]_

    metadata : dict
        The site metadata available in the file.

    See Also
    --------
    pvlib.iotools.parse_epw

    Notes
    -----

    The returned structures have the following fields.

    ===============   ======  =========================================
    key               format  description
    ===============   ======  =========================================
    loc               String  default identifier, not used
    city              String  site loccation
    state-prov        String  state, province or region (if available)
    country           String  site country code
    data_type         String  type of original data source
    WMO_code          String  WMO identifier
    latitude          Float   site latitude
    longitude         Float   site longitude
    TZ                Float   UTC offset
    altitude          Float   site elevation
    ===============   ======  =========================================


    +-------------------------------+-----------------------------------------+
    | EPWData field                 | description                             |
    +===============================+=========================================+
    | index                         | A pandas datetime index. NOTE, times are|
    |                               | set to local standard time (daylight    |
    |                               | savings is not included). Days run from |
    |                               | 0-23h to comply with PVLIB's convention.|
    +-------------------------------+-----------------------------------------+
    | year                          | Year, from original EPW file. Can be    |
    |                               | overwritten using coerce function.      |
    +-------------------------------+-----------------------------------------+
    | month                         | Month, from original EPW file.          |
    +-------------------------------+-----------------------------------------+
    | day                           | Day of the month, from original EPW     |
    |                               | file.                                   |
    +-------------------------------+-----------------------------------------+
    | hour                          | Hour of the day from original EPW file. |
    |                               | Note that EPW's convention of 1-24h is  |
    |                               | not taken over in the index dataframe   |
    |                               | used in PVLIB.                          |
    +-------------------------------+-----------------------------------------+
    | minute                        | Minute, from original EPW file. Not     |
    |                               | used.                                   |
    +-------------------------------+-----------------------------------------+
    | data_source_unct              | Data source and uncertainty flags. See  |
    |                               | [1]_, chapter 2.13                      |
    +-------------------------------+-----------------------------------------+
    | temp_air                      | Dry bulb temperature at the time        |
    |                               | indicated, deg C                        |
    +-------------------------------+-----------------------------------------+
    | temp_dew                      | Dew-point temperature at the time       |
    |                               | indicated, deg C                        |
    +-------------------------------+-----------------------------------------+
    | relative_humidity             | Relative humidity at the time indicated,|
    |                               | percent                                 |
    +-------------------------------+-----------------------------------------+
    | atmospheric_pressure          | Station pressure at the time indicated, |
    |                               | Pa                                      |
    +-------------------------------+-----------------------------------------+
    | etr                           | Extraterrestrial horizontal radiation   |
    |                               | recv'd during 60 minutes prior to       |
    |                               | timestamp, Wh/m^2                       |
    +-------------------------------+-----------------------------------------+
    | etrn                          | Extraterrestrial normal radiation recv'd|
    |                               | during 60 minutes prior to timestamp,   |
    |                               | Wh/m^2                                  |
    +-------------------------------+-----------------------------------------+
    | ghi_infrared                  | Horizontal infrared radiation recv'd    |
    |                               | during 60 minutes prior to timestamp,   |
    |                               | Wh/m^2                                  |
    +-------------------------------+-----------------------------------------+
    | ghi                           | Direct and diffuse horizontal radiation |
    |                               | recv'd during 60 minutes prior to       |
    |                               | timestamp, Wh/m^2                       |
    +-------------------------------+-----------------------------------------+
    | dir_n_rad                     | Amount of direct normal radiation       |
    |                               | (modeled) recv'd during 60 minutes prior|
    |                               | to timestamp, Wh/m^2                    |
    +-------------------------------+-----------------------------------------+
    | dif_h_rad                     | Amount of diffuse horizontal radiation  |
    |                               | recv'd during 60 minutes prior to       |
    |                               | timestamp, Wh/m^2                       |
    +-------------------------------+-----------------------------------------+
    | global_hor_illum              | Avg. total horizontal illuminance recv'd|
    |                               | during the 60 minutes prior to          |
    |                               | timestamp, lx                           |
    +-------------------------------+-----------------------------------------+
    | direct_normal_illum           | Avg. direct normal illuminance recv'd   |
    |                               | during the 60 minutes prior to          |
    |                               | timestamp, lx                           |
    +-------------------------------+-----------------------------------------+
    | diffuse_horizontal_illum      | Avg. horizontal diffuse illuminance     |
    |                               | recv'd during the 60 minutes prior to   |
    |                               | timestamp, lx                           |
    +-------------------------------+-----------------------------------------+
    | zenith_luminance              | Avg. luminance at the sky's zenith      |
    |                               | during the 60 minutes prior to          |
    |                               | timestamp, cd/m^2                       |
    +-------------------------------+-----------------------------------------+
    | wind_direction                | Wind direction at time indicated,       |
    |                               | degrees from north (360 = north; 0 =    |
    |                               | undefined,calm)                         |
    +-------------------------------+-----------------------------------------+
    | wind_speed                    | Wind speed at the time indicated, m/s   |
    +-------------------------------+-----------------------------------------+
    | total_sky_cover               | Amount of sky dome covered by clouds or |
    |                               | obscuring phenomena at time stamp,      |
    |                               | tenths of sky                           |
    +-------------------------------+-----------------------------------------+
    | opaque_sky_cover              | Amount of sky dome covered by clouds or |
    |                               | obscuring phenomena that prevent        |
    |                               | observing the sky at time stamp, tenths |
    |                               | of sky                                  |
    +-------------------------------+-----------------------------------------+
    | visibility                    | Horizontal visibility at the time       |
    |                               | indicated, km                           |
    +-------------------------------+-----------------------------------------+
    | ceiling_height                | Height of cloud base above local terrain|
    |                               | (7777=unlimited), meter                 |
    +-------------------------------+-----------------------------------------+
    | present_weather_observation   | Indicator for remaining fields: If 0,   |
    |                               | then the observed weather codes are     |
    |                               | taken from the following field. If 9,   |
    |                               | then missing weather is assumed.        |
    +-------------------------------+-----------------------------------------+
    | present_weather_codes         | Present weather code, see [1], chapter  |
    |                               | 2.9.1.28                                |
    +-------------------------------+-----------------------------------------+
    | precipitable_water            | Total precipitable water contained in a |
    |                               | column of unit cross section from earth |
    |                               | to top of atmosphere, cm. Note that some|
    |                               | old_TMY3.epw files may have incorrect   |
    |                               | unit if it was retrieved from           |
    |                               | www.energyplus.net.                     |
    +-------------------------------+-----------------------------------------+
    | aerosol_optical_depth         | The broadband aerosol optical depth per |
    |                               | unit of air mass due to extinction by   |
    |                               | aerosol component of atmosphere,        |
    |                               | unitless                                |
    +-------------------------------+-----------------------------------------+
    | snow_depth                    | Snow depth in centimeters on the day    |
    |                               | indicated, (999 = missing data)         |
    +-------------------------------+-----------------------------------------+
    | days_since_last_snowfall      | Number of days since last snowfall      |
    |                               | (maximum value of 88, where 88 = 88 or  |
    |                               | greater days; 99 = missing data)        |
    +-------------------------------+-----------------------------------------+
    | albedo                        | The ratio of reflected solar irradiance |
    |                               | to global horizontal irradiance,        |
    |                               | unitless                                |
    +-------------------------------+-----------------------------------------+
    | liquid_precipitation_depth    | The amount of liquid precipitation      |
    |                               | observed at indicated time for the      |
    |                               | period indicated in the liquid          |
    |                               | precipitation quantity field,           |
    |                               | millimeter                              |
    +-------------------------------+-----------------------------------------+
    | liquid_precipitation_quantity | The period of accumulation for the      |
    |                               | liquid precipitation depth field, hour  |
    +-------------------------------+-----------------------------------------+


    References
    ----------

    [1] `EnergyPlus documentation, Auxiliary Programs
       <https://energyplus.net/documentation>`
    [2] pvlib.iotools.parse_epw
    '''

    # Assume it's accessible via the file system
    csvdata = open(str(filename), 'r')
    try:
        data, meta = parse_epw(csvdata, coerce_year)
    finally:
        csvdata.close()
    return data, meta


def parse_epw(csvdata, coerce_year=None):
    """
    Given a file-like buffer with data in Energy Plus Weather (EPW) format,
    parse the data into a dataframe.

    Parameters
    ----------
    csvdata : file-like buffer
        a file-like buffer containing data in the EPW format

    coerce_year : None or int, default None
        If supplied, the year of the data will be set to this value. This can
        be a useful feature because EPW data is composed of data from
        different years.
        Warning: EPW files always have 365*24 = 8760 data rows;
        be careful with the use of leap years.

    Returns
    -------
    data : DataFrame
        A pandas dataframe with the columns described in the table
        below. For more detailed descriptions of each component, please
        consult the EnergyPlus Auxiliary Programs documentation
        available at: https://energyplus.net/documentation.

    metadata : dict
        The site metadata available in the file.

    See Also
    --------
    pvlib.iotools.read_epw
    """
    # Read line with metadata
    firstline = csvdata.readline()

    head = ['loc', 'city', 'state-prov', 'country', 'data_type', 'WMO_code',
            'latitude', 'longitude', 'TZ', 'altitude']
    meta = dict(zip(head, firstline.rstrip('\n').split(",")))

    meta['altitude'] = float(meta['altitude'])
    meta['latitude'] = float(meta['latitude'])
    meta['longitude'] = float(meta['longitude'])
    meta['TZ'] = float(meta['TZ'])

    colnames = ['year', 'month', 'day', 'hour', 'minute', 'data_source_unct',
                'temp_air', 'temp_dew', 'relative_humidity',
                'atmospheric_pressure', 'etr', 'etrn', 'ghi_infrared', 'ghi',
                'dir_n_rad', 'dif_h_rad', 'global_hor_illum',
                'direct_normal_illum', 'diffuse_horizontal_illum',
                'zenith_luminance',
                'wind_direction', 'wind_speed', 'total_sky_cover',
                'opaque_sky_cover', 'visibility', 'ceiling_height',
                'present_weather_observation', 'present_weather_codes',
                'precipitable_water', 'aerosol_optical_depth', 'snow_depth',
                'days_since_last_snowfall', 'albedo',
                'liquid_precipitation_depth', 'liquid_precipitation_quantity']

    # We only have to skip 6 rows instead of 7 because we have already used
    # the realine call above.
    data = pd.read_csv(csvdata, skiprows=6, header=0, names=colnames)

    # Change to single year if requested
    if coerce_year is not None:
        data["year"] = coerce_year

    # create index that supplies correct date and time zone information
    dts = data[['month', 'day']].astype(str).apply(lambda x: x.str.zfill(2))
    hrs = (data['hour'] - 1).astype(str).str.zfill(2)
    dtscat = data['year'].astype(str) + dts['month'] + dts['day'] + hrs
    idx = pd.to_datetime(dtscat, format='%Y%m%d%H')
    idx = idx.dt.tz_localize(int(meta['TZ'] * 3600))
    data.index = idx

    return data, meta


def sol_rad_tilt_surf(weather_data, surface_orientation, albedo):
    """
    Created on Fri Sep 10 11:04:48 2021
    @author: cghiaus

    Calculate solar radiation on a tilted surface from weathear data obtained
    from `*.epw` file.

    Parameters
    ----------
    weather_data : DataFrame
        Index : datetime64
        Column names :
            'temp_air' : dry bulb temperature at the time indicated, °C
            'dir_n_rad' : direct normal radiation during last 60 min, Wh/m²
            'dif_h_rad' : diffuse horizontal rad. during last 60 min, Wh/m²

    surface_orientation : dictionary
        'slope' : slope or tilt angle in deg: [0 180];
                    90°- vertical; > 90°- downward facing
        'azimuth' : surface azimuth in deg: [-180 180];
                    0-south; west-positive
        'latitude' : local latitude in deg: [-90 90],
                    north positive, south negative

    albedo : float
        diffuse reflection of solar radiation

    Returns
    -------
    solar_rad_tilt : DataFrame
        Index : datetime64
        Column names :
            'direct' : direct radiation on the surface, Wh/m²
            'diffuse' : diffuse radiation on the surface, Wh/m²
            'reflected' : reflected radiation on the surface, Wh/m²

    References
    ----------

    1. [Duffie 2020] J.A. Duffie, W. A. Beckman, N. Blair (2020) Solar
    Engineering of Thermal Processes, 5th ed. John Wiley & Sons, Inc.
    ISBN 9781119540281

    2. [Th-CE 2005] Réglementation Thermique 2005. Méthode de calcul Th-CE.
    Annexe à l’arrêté du 19 juillet 2006
    """
    B = surface_orientation['slope']
    Z = surface_orientation['azimuth']
    L = surface_orientation['latitude']

    # Transform degrees in radians
    B = B * np.pi / 180
    Z = Z * np.pi / 180
    L = L * np.pi / 180

    n = weather_data.index.dayofyear

    # [Duffie 2020] eq. 1.6.1a
    # [Th-CE] §11.2.1.1, eq. (78)
    declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
    d = declination_angle * np.pi / 180

    # [Duffie 2020] Example 1.6.1
    hour = weather_data.index.hour
    minute = weather_data.index.minute + 60
    hour_angle = 15 * ((hour + minute / 60) - 12)
    h = hour_angle * np.pi / 180

    # [Duffie 2020] incidence angle eq. 1.6.2
    # [Th-CE 2005] §11.2.1.1
    theta = np.sin(d) * np.sin(L) * np.cos(B)
    theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
    theta += np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
    theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
    theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
    theta = np.array(np.arccos(theta))
    theta[theta > (np.pi / 2)] = np.pi / 2

    # Direct radiation on a wall
    # [Th-CE 2005] §11.2.1.1
    dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
    dir_rad[dir_rad < 0] = 0

    # Diffuse radiation on a wall
    # [Th-CE 2005] §11.2.1.2, Eq. 79, p. 31
    dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(B)) / 2

    # Solar radiation reflected by the ground
    # [Th-CE 2005] §112.1.3, after eq. (78)
    gamma = np.cos(d) * np.cos(L) * np.cos(h)
    gamma += np.sin(d) * np.sin(L)
    gamma = np.array(np.arcsin(gamma))
    gamma[gamma < 1e-5] = 1e-5

    # Radiation reflected by the ground
    # [Th-CE 2005] §11.2.1.3 eq. (80)
    # direct radiation on horizontal surface
    dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)
    # total reflected radiation
    ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
    # reflected radiation eceived by the wall
    ref_rad *= (1 - np.cos(B) / 2)

    solar_rad_tilt = pd.DataFrame({'direct': dir_rad,
                                   'diffuse': dif_rad,
                                   'reflected': ref_rad})
    return solar_rad_tilt


def wall2TC(wall_types, walls_def, prefix="w"):
    """
    Created on Tue Jun 13 17:19:23 2023
    based on wall2TC()
    @author: cghiaus

    Creates a thermal circuit TC as a set of differential-algebraic-equations
    (DAE) A, G, C, b, f, y for internal walls
    from the files `wall_type` and `walls_def`.

    Parameters
    ----------
    wall_types : str
        Name of .csv file describing the types of the walls
    walls_def : str
        Name of .csv file giving the data on walls.
    prefix : str
        Prefix of the ID of the thermal circuit in indexes and columns.

    Returns
    -------
    TC : dict
        Thermal circuit: A, G, C, b, f, y.
        The index k of TC[k] represents the wall ID in walls_data.csv file

    Descrition of input files
    -------------------------
    wall_types.csv: example
        type, Material, Conductivity, Specific heat, Density, Width, Mesh\n
        0, Concrete,     1.4,  880, 2300,  0.2, 2\n
        0, Insulation, 0.027, 1210,   55, 0.08, 1\n
        1, Glass,        1.4,  750, 2500,0.004, 1\n

        type : str
            type of the wall; same type on more rows for multilayer wall
        Material : str
            material
        Conductivity : float
            thermal conductivity of material, W·m⁻¹·K⁻¹
        Specific heat : float
            specific heat capacity of material, J⋅kg⁻¹⋅K⁻¹
        Density : float
            kg⋅m⁻³
        Width : float
            m
        Mesh : int
            number of meshes in space discretization of material

    walls_def.csv:
        3 kinds of walls:
            - generic: T0, T1 specified or not
            - out: T0 specified, without T1
            - in: without T0 and T1
        ID,type,A,β,Q0,Q1,h0,h1,α0,α1,ε0,ε1,y\n
        0,0,10,90,0,Qo,Qi,25,8,0.25,0.30,0.85,0.70,"[0, -1]"\n

        ID : str
            wall instance identifiant
        type : str
            value from wall_types.csv
        Area : float
            surface area of the plan wall, m²
        [β] : float or nan; for generic and out walls
            wall slope, °; 90° vertical; > 90° downward
        [γ] : float or nan; for generic and out walls
            azimuth, °; 0° South, 180° North, >0 westward, <0 eastward
        [T0] : str or nan; for generic and out walls
            name of the temperature source of surface 0, °C
        [T1] : str or nan; only for generic walls
            name of the temperature source of surface 1, °C
        Q0 : str or NaN
            name of the flow rate source of surface 0, W
        Q1 : str or nan
            name of the flow rate source of surface 1, W
        h0 : float
            convection coefficient surface 0, W·m⁻²·K⁻¹
        h1 : float
            convection coefficient surface 1, W·m⁻²·K⁻¹
        α0 : float [0, 1] or nan
            short-wave absorbtion coefficient of surface 0, -
        α1 : float [0, 1] or nan
            short-wave absorbtion coefficient of surface 1, -
        ε0 : float [0, 1] or nan
            long-wave hemispherical emissivity of surface 0, -
        ε1 : float [0, 1] or nan
            long-wave hemispherical emissivity of surface 1, -
        y : slicing indexer or nan
            output temperature nodes by using slicing, e.g. [0, -1]


    Description of outputs
    ----------------------
    TC[k]:
        TC[k]['A']: Dataframe
            incidence matrix A of wall k
        TC[k]['G']: Series
            conductance vector G of wall k
        TC[k]['C']: Series
            capacity vector C of wall k
        TC[k]['b']: Series
            temperature source vector b of wall k
        TC[k]['f']: Series
            flow rate source vector f of wall k
        TC[k]['y']: Series
            output temperature nodes vector: 1 for output, 0 otherwise

    In TC[k], the indexes are:
        - {name}. the wall {name} followed by number k
        - q.. the flow number of wall k
        - θ.. the temperature number of wall k
        Example: w0_q0 for flows (branches) or w0_θ0 for temperatures (nodes)
    """

    def wall_instances(wall_types, walls_def):
        """
        Merges `walls_def` with `wall_types` to obtain a df for all layers
        of the walls.

        Parameters
        ----------
        wall_types : DataFrame
            Composition of walls and physical properties for
            materials of each layer.
        walls_def : DataFrame
            Definition of each wall instance. Three kind of wall definitions:
            generic, outdoor & indoor walls

        Returns
        -------
        wall : DataFrame
            Tidy DataFrame walls and layers in walls.
        """
        wall = walls_def.merge(wall_types, on='type')
        wall = wall.sort_values("ID")
        # conductances
        wall['U'] = wall['Conductivity'] / wall['Width']
        # capacities
        volume = wall['Width'] * wall['Area']
        wall['C'] = wall['Density'] * wall['Specific heat'] * volume
        return wall

    def number_branches(wall):
        """
        Gives the number of branches in each wall as a function of number of
        layers and number of meshes in each layer of a wall. Each mesh has
        two resistances (corresponding to 1/2 of the width of the mesh) and
        one capacity (corresponding to the whole volume of the mesh).

        Parameters
        ----------
        wall : DataFrame
            Tidy DataFrame of walls and layers in walls.

        Returns
        -------
        nq : Series
            Number of branches for each wall in walls.

        """
        # 1 mesh: 2 R and 1 C
        mesh_df = wall.groupby("ID").agg({"Mesh": "sum"})
        mesh = mesh_df.squeeze()   # Series or int <-- DataFrame
        if not isinstance(mesh, pd.Series):
            mesh = mesh_df['Mesh'].rename_axis(None)
        # number of flow branches in each wall:
        # 2R per mesh for conduction & 2 convection for each wall face
        nq = (2 * (mesh + 1)).rename('nq')
        return nq

    def DAE_without_bound_temp(walls_def, wall, nq, k):
        """
        Creates the matrices A, G, C and vectors b, f of the sysem of
        Differential Algebraic Equations (DAE) by not considering the
        temperature sources T0 and T1 (on the boundaries).

        Parameters
        ----------
        walls_def : DataFrame
            Definition of each wall instance. Three kind of wall definitions:
            generic, outdoor & indoor walls.
        wall : DataFrame
            Tidy DataFrame of walls and layers in walls.
        nq : Series
            Number of branches for each wall in walls.
        k : str
            ID of the wall from walls_def.

        Returns
        -------
        A, C : Array
        G, b, f : Series
        """
        # Boundary conditions: [To, Ti] = [NaN, NaN]
        A = np.diff(np.eye(nq[k] + 1)).T

        # G
        # conduction: #R = 2 * #mesh
        Um = wall.loc[wall["ID"] == k][['U', 'Mesh']]   # λ/w for mesh

        # mesh U: U -> [U/(2*mesh), U/(2*mesh), ...]
        U = Um['U'] * (2 * Um['Mesh'])  # value * (2*mesh) -> series
        repeat_index = np.repeat(Um.index, 2 * Um.Mesh)
        U = U.reindex(repeat_index)

        # append convection
        h0 = walls_def.loc[walls_def['ID'] == k, 'h0'].values[0]
        h1 = walls_def.loc[walls_def['ID'] == k, 'h1'].values[0]
        U = pd.concat([pd.Series([h0]),
                       U,
                       pd.Series([h1])], ignore_index=True)
        area = walls_def.loc[walls_def['ID'] == k, 'Area'].values[0]
        G = U * area

        # C
        Cm = wall.loc[wall["ID"] == k][['C', 'Mesh']]   # ρ·c·w·A for mesh

        # mesh C: G -> [C/mesh, ...]
        C = np.zeros([1, 2 * sum(Cm.Mesh) + 1])
        Cv = Cm.C.div(Cm.Mesh)            # value / mesh -> series
        repeat_index = np.repeat(Cm.index, Cm.Mesh)
        Cv = Cv.reindex(repeat_index)

        # insert zeros
        C = np.array(list(zip(np.zeros(Cv.shape[0]), Cv))).flatten()
        C = np.append(C, 0)     # last zero

        # insert zeros for Boundary Conditions: [nan, nan]
        C = np.insert(C, obj=0, values=0)
        C = np.append(C, 0)

        # b
        b = np.zeros([A.shape[0]])
        b = pd.Series(b)

        # f
        f = np.zeros([A.shape[1]])
        f = pd.Series(f)
        f.iloc[1] = walls_def.loc[walls_def['ID'] == k, 'Q0'].values[0]
        f.iloc[-2] = walls_def.loc[walls_def['ID'] == k, 'Q1'].values[0]

        return A, G, C, b, f

    def DAE_with_bound_temp(walls_def, A, G, C, b, f, k):
        """
        Add boundary conditions:
            - [To, nan]
            in A, C, f, y; delete first column: 0
            in b: insert To
            - [nan, Ti]
            in A, C, f; delete last column: -1
            in b: insert Ti
            - [To, Ti]
            in A, C, f; delete first and last column: [0, -1]
            in b: insert [To, Ti]

        Parameters
        ----------
        walls_def : DataFrame
            Definition of each wall instance.

        A, C : Array

        G, b, f : Series

        k : str
            ID of the wall from walls_def.
        """
        if 'T0' in walls_def.columns and 'T1' in walls_def.columns:
            """
            General wall with columns for sources T0 and T1
            """
            # Add temperature boundary conditions
            bc_T0 = walls_def.loc[walls_def['ID'] == k, ['T0']].notna().values
            bc_T1 = walls_def.loc[walls_def['ID'] == k, ['T1']].notna().values

            if bc_T0 and (not(bc_T1)):
                A = np.delete(A, 0, axis=1)
                C = np.delete(C, 0)
                b.iloc[0] = walls_def.loc[
                    walls_def['ID'] == k, 'T0'].values[0]
                f = f[1:]

            if not(bc_T0) and bc_T1:
                A = np.delete(A, -1, axis=1)
                C = np.delete(C, -1)
                b.iloc[-1] = '-' + walls_def.loc[
                    walls_def['ID'] == k, 'T1'].values[0]
                f = f[:-1]

            if bc_T0 and bc_T1:
                A = np.delete(A, [0, -1], axis=1)
                C = np.delete(C, [0, -1])
                b.iloc[0] = walls_def.loc[
                    walls_def['ID'] == k, 'T0'].values[0]
                # -Ti since the flow is from + to - in source Ti
                b.iloc[-1] = '-' + walls_def.loc[
                    walls_def['ID'] == k, 'T1'].values[0]
                f = f[1:-1]

        elif 'T0' in walls_def.columns and 'T1' not in walls_def.columns:
            """
            Outdoor wall with T0 (source out)
            """
            A = np.delete(A, 0, axis=1)
            C = np.delete(C, 0)
            b.iloc[0] = walls_def.loc[
                walls_def['ID'] == k, 'T0'].values[0]
            f = f[1:]

        return A, G, C, b, f

    def DAE_output(walls_def, A):
        """
        Add output vector `y`. The elements of `y` are `1`if the temperature
        node is an output and `0` otherwise.

        Parameters
        ----------
        walls_def : DataFrame
            Definition of each wall instance.
        A : Array
            Incidence matrix.

        Returns
        -------
        y : Series
            1 if the node is an output, 0 otherwise.

        """
        y = np.zeros([A.shape[1]])
        y = pd.Series(y)

        slice_str = walls_def.loc[walls_def['ID'] == k, 'y'].values[0]
        if type(slice_str) == str:
            parsed_slice = ast.literal_eval(slice_str)
            y.iloc[parsed_slice] = 1
        return y

    def DAE_pd(A, G, C, b, f, y, k):
        """
        Converts A, G, C, b, f, y of the DAE model into DataFrame (for A) and
        Series (for G, C, b, f, y) indexed with:
            - the ID of the wall from `wall_def`;
            - `θ`for temperature nodes and `q` for flow branches;
            - number of node of branch.
        e.g., `w3_θ1` , `w3_q1`

        Parameters
        ----------
        A, C : Array

        G, b, f, y : Series

        k : str
            ID of the wall from walls_def.

        Returns
        -------
        A : DataFrame

        G, C, b, f, y : Series
        """
        # DataFrame with wall-ID, q-flows & θ-tempertaures
        w_q = [prefix + str(k) + '_q' + str(x) for x in range(A.shape[0])]
        w_θ = [prefix + str(k) + '_θ' + str(x) for x in range(A.shape[1])]

        A = pd.DataFrame(data=A,
                         index=w_q,
                         columns=w_θ)
        G = pd.Series(G)
        C = pd.Series(C)

        G = G.set_axis(w_q)
        C = C.set_axis(w_θ)
        b = b.set_axis(w_q)
        f = f.set_axis(w_θ)
        y = y.set_axis(w_θ)

        return A, G, C, b, f, y

    walls = wall_instances(wall_types, walls_def)
    nq = number_branches(walls)

    TC = {}
    for k in walls_def['ID']:
        A, G, C, b, f = DAE_without_bound_temp(walls_def, walls, nq, k)
        A, G, C, b, f = DAE_with_bound_temp(walls_def, A, G, C, b, f, k)
        y = DAE_output(walls_def, A)
        A, G, C, b, f, y = DAE_pd(A, G, C, b, f, y, k)

        TC[prefix + k] = {'A': A,
                          'G': G,
                          'C': C,
                          'b': b,
                          'f': f,
                          'y': y}
    return TC


def file2TC(TC_file, name="w", auto_number=False):
    """
    Created on Wed Nov 16 12:02:41 2022
    @author: cghiaus

    Creates a thermal circuit TC as a set of differential-algebraic-equations
    (DAE) A, G, C, b, f from the file `TC_file`.

    Parameters
    ----------
    TC_file : str
        Name of .csv file describing the types of the walls
    name : str
        Name of the thermal circuit that will appear in the indexes of
        the returned thermal circuit TC.
    auto_number : bool
        Automatic numbering of indexes and columns

    Returns
    -------
    TC : dict
    TC = {"A": DataFrame,
          "G": Series,
          "C": Series,
          "b": Series,
          "f": Series
          "y": Series}
    Thermal circuit: A, G, C, b, f, y.

    A indexes: name_θ0, name_θ1, ..., columns: name_q0, name_q1, ...
    G indexes: name_q0, name_q1, ...
    C indexes: name_θ0, name_θ1, ...
    b indexes: name_q0, name_q1, ...
    f indexes: name_θ0, name_θ1, ...
    y indexes: name_θ0, name_θ1, ...

    Description of `TC_file`
    -----------------------
    TC_file.csv (file with `NaN`):
        A	θ0	θ1	θ2	θ3	θ4	G	b
        q0	1					250	To
        q1	-1	1				140
        q2		-1	1			140
        q3			-1	1		6.75
        q4				-1	1	6.75
        q5					-1	80	-Ti
        C		4048000		53240
        f	Qo				Qi

    or (file with `0`):
        A	θ0	θ1	θ2	θ3	θ4	  G	b
        q0	1	0	0	0	0	  250	To
        q1	-1	1	0	0	0	  140	0
        q2	0	-1	1	0	0	  140	0
        q3	0	0	-1	1	6.75    0  	0
        q4	0	0	0	-1	1	  6.75	0
        q5	0	0	0	0	-1	   80	-Ti
        C	0	4048000	0	53240	0	0	0
        f	Qo	0	0	0	Qi	    0	0
    """
    TC_file = pd.read_csv(TC_file, index_col=0)
    TC_file.index.name = None
    TC_file = TC_file.fillna(0)

    # select A, G, C, b, f, y from TC_file
    A = TC_file.iloc[:-3, :-2].astype(float)
    G = TC_file.iloc[:-3, -2].astype(float)
    C = TC_file.iloc[-3, :-2].astype(float)
    b = TC_file.iloc[:-3, -1]
    f = TC_file.iloc[-2, :-2]
    y = TC_file.iloc[-1, :-2].astype(int)

    if auto_number:
        """
        Auto-numbering of the indexes and columns with _q and _θ.
        """
        # Indexes wall-number, q-flows & θ-temperatures
        tc_q = [name + '_q' + str(x) for x in range(A.shape[0])]
        tc_θ = [name + '_θ' + str(x) for x in range(A.shape[1])]
    else:
        tc_q = [name + '_' + str(x) for x in list(A.index)]
        tc_θ = [name + '_' + str(x) for x in list(A.columns)]

    A.index, A.columns = tc_q, tc_θ
    G.index = tc_q
    C.index = tc_θ
    b.index = tc_q
    f.index = tc_θ
    y.index = tc_θ

    # Thermal circuit
    TC = {"A": A,
          "G": G,
          "C": C,
          "b": b,
          "f": f,
          "y": y}
    return TC


def bldg2TCd(folder_path, TC_auto_number):
    """
    Created on Tue Jun 20 19:56:31 2023

    @author: cghiaus

    Convert a folder containing the files characterizing the building into
    a disassambled thermal circuit.

    Calls :

    wall2TC() : from wall_types.csv and (wall_in.csv, wall_out.cvs,
    wall_generic.cvs) to thermal circuit TC = (A, G, C, b, f, y).

    file2TC() : from file TC_.csv to thermal circuit TC = (A, G, C, b, f, y).

    Parameters
    ----------
    bldg_path : str
        Path of the folder containing *.csv files describing the disassambled
        thermal cicuit: walls, thermal circuits, and assembly matrix and/or
        lists.

    Returns
    -------
    TCd : dict
        Disassembled thermal circuits. The number of the thermal circuit TC
        in the disassembled thermal circuit TCd is the number of walls in file
        `walls_in.csv`, `walls_out.csv`, `walls_generic.csv` plus the number of
        files describing the thermal circuits TC_.csv.

        Indexes of TC in TCd:
            - w for walls; e.g. w2_n1 for wall 2 node 1
            - c for thermal circuits; e.g. c1_b2 for circuit TC1 branch 2

        Each circuit is a dictionary:
            - A: DataFrame, incidence matrix;
            - G: Series, diagonal of the conductance matrix;
            - C: Series, diagonal of the capacity matrix;
            - b: Series, vector of temperature sources on branches;
            - f: Series, vector of flow sources in nodes;
        with indexes `b` for branches, e.g. w0_b1 for wall_0 branch_1,
        and `n` for nodes, e.g. c1_n2 for circuit TC1 node 2.

    Description of the folder containing the disassembled thermal circuit:
        - Assembly_matrix: pairs of nodes that are in common in the assembled
        circuit:
            - ID0 node0  ID1  node1
            - 0   4       1   0
            means node 4 of circuit 0 is in common with node 0 of circuit 1.
            Node 0 of circuit 1 is deleted after assembling.
        - Assembly_lists: [[[circuit, node], [circuit, node], [circuit, node]],
                           ...,
                           [[circuit, node], [circuit, node], [circuit, node]]
                           ]
        lists of lists of pairs [circuit, node] which are put in common

        - TC*.csv: thermal circui: see file2TC()

        - wall_types: see wall2TC()
        - walls_in, walls_out, walls_generic: see wall2TC()

    Description of indexes of thermal circuits TC:
        - walls: w#_, where # is the wall number (ID); e.g. w0_b0 or w0_n0
        - TC: c#_, where # is the TC number; e.g. c1_n0

    Indexes of branches and nodes:
        - branch indexes: b; e.g. w0_b0 for wall 0 branch 0
        - node indexes: n; e.g. c1_n2 for TC1 node 2.

    How to access of disassambled thermal circuit TCd:
        - a circuit: TCd[1]
        - an element of the circuit: TCd[1]['A']
    """

    file_path = os.path.join(folder_path, "wall_types.csv")
    wall_types = pd.read_csv(file_path)

    TCd = {}

    # Thermal circuits from walls_.csv & wall_types.csv files of walls
    file_path = os.path.join(folder_path, "walls_generic.csv")
    if os.path.isfile(file_path):
        walls = pd.read_csv(file_path)
        TCd_generic = wall2TC(wall_types, walls, prefix="g")
        TCd.update(TCd_generic)

    file_path = os.path.join(folder_path, "walls_in.csv")
    if os.path.isfile(file_path):
        walls = pd.read_csv(file_path)
        TCd_in = wall2TC(wall_types, walls, prefix="i")
        TCd.update(TCd_in)

    file_path = os.path.join(folder_path, "walls_out.csv")
    if os.path.isfile(file_path):
        walls = pd.read_csv(file_path)
        TCd_out = wall2TC(wall_types, walls, prefix="o")
        TCd.update(TCd_out)

    # Read all files from folder with data on building
    csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

    TC_files = [TC_file for TC_file in csv_files if 'TC' in TC_file]
    TC_files = sorted(TC_files)

    # Thermal circuits from files of thermal circuit
    for k in range(len(TC_files)):
        name = "c" + str(k)
        TCd[name] = file2TC(TC_files[k], name, TC_auto_number)

    return TCd


# def assemble_TCd(folder_path):
#     # Obtain dissambled TCd from files in folder_path
#     TCd = bldg2TCd(folder_path)

#     # Assembly matrix from Assembly_matrix.csv in folder_path
#     csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
#     Assembly_file = [Asf for Asf in csv_files if 'Assembly_int' in Asf][0]
#     Ass = pd.read_csv(Assembly_file)

#     # Rename common columns in A, C, f (rows remain the same)
#     for k in range(len(Ass)):
#         old_name = TCd[Ass['ID1'][k]]['A'].columns[Ass['node1'][k]]
#         new_name = TCd[Ass['ID0'][k]]['A'].columns[Ass['node0'][k]]

#         TCd[Ass['ID1'][k]]['A'] = TCd[
#             Ass['ID1'][k]]['A'].rename(columns={old_name: new_name})

#         TCd[Ass['ID1'][k]]['C'] = TCd[
#             Ass['ID1'][k]]['C'].rename(index={old_name: new_name})

#         TCd[Ass['ID1'][k]]['f'] = TCd[
#             Ass['ID1'][k]]['f'].rename(index={old_name: new_name})

#     TC = {'A': pd.DataFrame(),
#           'G': pd.Series(dtype=float),
#           'C': pd.Series(dtype=float),
#           'b': pd.Series(dtype=str),
#           'f': pd.Series(dtype=str)}

#     for k in range(len(TCd)):
#         TC['A'] = TC['A'].add(TCd[k]['A'], fill_value=0).fillna(0)
#         TC['G'] = TC['G'].add(TCd[k]['G'], fill_value=0)
#         TC['C'] = TC['C'].add(TCd[k]['C'], fill_value=0)

#         # b and f are Series of str
#         TCd[k]['b'] = TCd[k]['b'].replace(0, '')
#         TC['b'] = pd.concat([TC['b'], TCd[k]['b']])
#         TCd[k]['f'] = TCd[k]['f'].replace(0, '')
#         TC['f'] = pd.concat([TC['f'], TCd[k]['f']])

#     # sum the same indexes in TC['f']
#     TC['f'] = TC['f'].groupby(TC['f'].index).sum()

#     # Sort
#     TC['A'] = TC['A'].sort_index()
#     TC['G'] = TC['G'].sort_index()
#     TC['C'] = TC['C'].sort_index()
#     TC['b'] = TC['b'].sort_index()
#     TC['f'] = TC['f'].sort_index()

#     return TC


def assemble_TCd_matrix(TCd, ass_mat, TC_auto_number=False):
    """

    Parameters
    ----------
    TCd : dict
        Dictionary of disassembled thermal circuits: A, G, C, b, f, y. Each
        circuit has an ID.

    ass_mat : DataFrame
        Assembly matrix. Four columns: TC0, node 0, TC1, node 1
        of the nodes that are assembled. The nodes are described by the ID
        of the thermal circuit, TC0 and TC1, and the number of the node,
        respecting the convention of indexes for list, e.g., 0 for 1st element,
        -1 for the last element.

    Returns
    -------
    TC : dict
        Thermal circuit: A, G, C, b, f, y.

    To obtain ass_mat, read the file:
    assembly_matrix.csv

    TC0,    node0,  TC1,    node1
    ow0,    4,      c0,     0
    c1,     -2,     c0,     1
    c2,     0,      ow0,    -1
    c2,     0,      c3,     0
    c2,     0,      c1,     -1
    """

    def find_index_name(TC_name, row_index):
        TC = TCd.get(TC_name)
        if TC:
            series = TC.get('C')
            if series is not None:
                if row_index >= 0:
                    return series.index[row_index]
                else:
                    return series.index[row_index % len(series)]
        return None

    # Rename common columns in A, C, f (rows remain the same)
    for index, row in ass_mat.iterrows():
        name_old = find_index_name(TC_name=row['TC1'], row_index=row['node1'])
        name_new = find_index_name(TC_name=row['TC0'], row_index=row['node0'])
        # print(name_new, "<--", name_old)

        TCd[row['TC1']]['A'] = TCd[
            row['TC1']]['A'].rename(columns={name_old: name_new})
        TCd[row['TC1']]['C'] = TCd[
            row['TC1']]['C'].rename(index={name_old: name_new})
        TCd[row['TC1']]['f'] = TCd[
            row['TC1']]['f'].rename(index={name_old: name_new})
        TCd[row['TC1']]['y'] = TCd[
            row['TC1']]['y'].rename(index={name_old: name_new})

    # Add the nodes with the same name
    TC = {'A': pd.DataFrame(),
          'G': pd.Series(dtype=float),
          'C': pd.Series(dtype=float),
          'b': pd.Series(dtype=str),
          'f': pd.Series(dtype=str),
          'y': pd.Series(dtype=int)}

    for k in TCd.keys():
        TC['A'] = TC['A'].add(TCd[k]['A'], fill_value=0).fillna(0)
        TC['G'] = TC['G'].add(TCd[k]['G'], fill_value=0)
        TC['C'] = TC['C'].add(TCd[k]['C'], fill_value=0)

        # b and f are Series of str
        TCd[k]['b'] = TCd[k]['b'].replace(0, '')
        TC['b'] = pd.concat([TC['b'], TCd[k]['b']])
        TCd[k]['f'] = TCd[k]['f'].replace(0, '')
        TC['f'] = pd.concat([TC['f'], TCd[k]['f']])

        # y
        TC['y'] = TC['y'].add(TCd[k]['y'], fill_value=0)

    # sum the same indexes in TC['f']
    TC['f'] = TC['f'].groupby(TC['f'].index).sum()

    # Sort
    TC['A'] = TC['A'].sort_index()
    TC['G'] = TC['G'].sort_index()
    TC['C'] = TC['C'].sort_index()
    TC['b'] = TC['b'].sort_index().replace('', 0)
    TC['f'] = TC['f'].sort_index().replace('', 0)
    TC['y'] = TC['y'].sort_index().replace('', 0)

    return TC


def assemble_lists2matrix(ass_lists):
    """
    Converts ass_lists into ass_matrix (to be used by assemble_TCd_matrix()).

    Parameters
    ----------
    ass_list : TYPE
        DESCRIPTION.

    Returns
    -------
    ass_mat : DataFrame
        Assembly matrix. Four columns: TC0, node 0, TC1, node 1
        of the nodes that are assembled. The nodes are described by the ID
        of the thermal circuit, TC0 and TC1, and the number of the node,
        respecting the convention of indexes for list, e.g., 0 for 1st element,
        -1 for the last element.

    For example:
    To obtain ass_lists DataFrame, read the file:

    ass_lists.csv
    node0,          nodes
    "['c2', 0]",    "['ow0', -1], ['c3', 0], ['c1', -1]"
    "['ow0', 4]",   "['c0', 0],"
    "['c1', -2]",   "['c0', 1],"

    with:
    ass_lists = pd.read_csv(folder_path + '/assembly_lists.csv')

    Note: Python accepts the character ' of " for str.
    """

    # Initialize an empty list to store the transformed data
    ass_mat_data = []

    # Iterate over each row in the DataFrame
    for index, row in ass_lists.iterrows():
        # Convert the string representation of the list to an actual list
        node0 = ast.literal_eval(row['node0'])

        # Iterate over each list in the nodes column
        for sublist in ast.literal_eval(row['nodes']):
            # Extract the values from the sublist
            tc1, node1 = sublist

            # Append the transformed data to the list
            ass_mat_data.append(
                {'TC0': node0[0],
                 'node0': node0[1],
                 'TC1': tc1,
                 'node1': node1})

    # Create the transformed DataFrame
    ass_mat = pd.DataFrame(ass_mat_data)

    return ass_mat


def tc2ss(TC):
    """
    Implements alogithm from
    Ghiaus, C. (2013). Causality issue in the heat balance method for
    calculating the design heating and cooling load. Energy, 50, 292-301.

    Parameters
    ----------
    TC : dict
        keys: 'A', 'C', 'G', 'b', 'f', 'y'
    TC['A'] : DataFrame
        arc-node incidence matrix:
        index: heat flow rates
        cols: temperature nodes
        data: -1 if flow gets out; 1 if flow gets in; 0 if flow not in/out

    TC['G'] : Series
        diagonal of matrix of conductances
        index: heat flow rates
        data: conductances or 0

    TC['C'] : Series
        diagonal of matrix of capacities
        index: temperature nodes
        data: capacities or 0

    TC['b'] : Series
        name of temperature sources on branches
        index: temperature nodes
        data:  temperature source name (TYPE str) or '' (empty string)

    TC['f'] : Series
        name of flow-rate sources in nodes
        index: heat flow rates
        data: flow-rate source name (TYPE str) or '' (empty string)

    TC['y'] : Series of float64
        vector indicating the ouput nodes
        index: temperature nodes
        data: 0 if not an output; int > 0 if output

    Returns
    -------
    As state matrix in state equation
    Bs input matrix in state equation
    Cs output matrix in observation equation
    Ds input matrix in observation equation
    us list of inputs
    """

    def inv(A):
        """
        Inverse of a matrix contained in DataFrame.
        The index and column names are kept.

        Parameters
        ----------
        A : DataFrame
            a matrix.

        Returns
        -------
        A_inv : DataFrame
            inverse of matrix.

        """
        A_inv = pd.DataFrame(np.linalg.inv(A.values),
                             index=A.index,
                             columns=A.columns)
        return A_inv

    A, G, C, b, f, y = TC['A'], TC['G'], TC['C'], TC['b'], TC['f'], TC['y']
    
    # indexes of inputs
    idx_G = pd.Series(G.index)
    idx_C = pd.Series(C.index)
    idx_u = idx_G.append(idx_C)
    
    G = pd.DataFrame(np.diag(G), index=G.index, columns=G.index)
    C = pd.DataFrame(np.diag(C), index=C.index, columns=C.index)

    C_diag = np.diag(C.values)  # Get the values on the diagonal of C

    # Input vector, Eq.(22)
    f0 = f.loc[C_diag == 0]
    fc = f.loc[C_diag != 0]
    u = pd.concat([b, f0, fc])

    # DAE (Differential Algebraic Equations), Eq.(13)
    K = -A.T @ G @ A
    Kb = A.T @ G

    # Partition K and Kb based on values of diagonal of C, Eq.(14)
    K00 = K.loc[C_diag == 0, C_diag == 0]
    K01 = K.loc[C_diag == 0, C_diag != 0]
    K10 = K.loc[C_diag != 0, C_diag == 0]
    K11 = K.loc[C_diag != 0, C_diag != 0]

    Kb0 = Kb.loc[C_diag == 0]
    Kb1 = Kb.loc[C_diag != 0]

    C0 = C.loc[C_diag == 0, C_diag == 0]
    Cc = C.loc[C_diag != 0, C_diag != 0]

    eye00 = pd.DataFrame(np.eye(len(C0)),
                         index=C0.index, columns=C0.columns)
    eye11 = pd.DataFrame(np.eye(len(Cc)),
                         index=Cc.index, columns=Cc.columns)
    zeros = pd.DataFrame(np.zeros((len(C0), len(Cc))),
                         index=C0.index, columns=Cc.columns)
    # State equation, Eqs.(20), (21)
    As = inv(Cc) @ (-K10 @ inv(K00) @ K01 + K11)
    Bs = inv(Cc) @ pd.concat([-K10 @ inv(K00) @ Kb0 + Kb1,
                              -K10 @ inv(K00),
                              eye11], axis=1)

    Bs = Bs.reindex(columns=idx_u)

    # Observation equation for all temperature nodes
    # - for nodes with capacities θc
    Csc = pd.DataFrame(np.eye(len(Cc)),
                       index=Cc.index, columns=Cc.columns)
    Dsc = pd.DataFrame(np.zeros((len(Cc), len(u))),
                       index=Cc.index, columns=u.index)
    # - for nodes without capacities θ0, Eq.(24), (26), (27)
    Cs0 = -inv(K00) @ K01
    Ds0 = -inv(K00) @ pd.concat([Kb0, eye00, zeros], axis=1)
    # -for all nodes: without capacity θ0, with capacity θc
    Cs = pd.concat([Cs0, Csc], axis=0)
    Ds = pd.concat([Ds0, Dsc], axis=0)
    Ds = Ds.reindex(columns=idx_u)

    u = u.reindex(idx_u)

    # From u, extract actual inputs (inputs that are not identically zero)
    non_zero_inputs = pd.concat([(TC['b'] != 0), (TC['f'] != 0)])
    Bs = Bs.loc[:, non_zero_inputs]
    Ds = Ds.loc[:, non_zero_inputs]
    u = u.loc[non_zero_inputs]

    # From complete output vector [θ0, θc], extract actual outputs given by y
    Cs = Cs.loc[y.loc[y != 0].index]
    Ds = Ds.loc[y.loc[y != 0].index]

    return As, Bs, Cs, Ds, u


def inputs_in_time(us, input_data):
    """
    Input vectors in time.

    Parameters
    ----------
    us : Series
        Index: name of the branch or the node.
        Value: 'str' name of the variable that needs to be in `data`
    data :
        Variables containing the inputs in time.

    Returns
    -------
    u : DataFrame
        rows: inputs values at a moment of time
        columns: variables from `data`in the order from `us`.

    Example:
        given
    data = {'To': To, 'Ti_sp': Ti_sp, 'Φ': Φ, 'Qa': Qa, 'Qo': Qo, 'Qi': Qi}
        call

    u = inputs_in_time(us, input_data)
    """
    input_data = pd.DataFrame(input_data)
    u = input_data[us.values]
    return u


def print_TC(TC):
    """

    Parameters
    ----------
    TC : dict
        Thermal circuit.

    Returns
    -------
    None.

    """
    print('A:')
    print(TC['A'], '\n')
    print('G:')
    print(TC['G'], '\n')
    print('C:')
    print(TC['C'], '\n')
    print('b:')
    print(TC['b'], '\n')
    print('f:')
    print(TC['f'], '\n')
    print('y:')
    print(TC['y'], '\n')
