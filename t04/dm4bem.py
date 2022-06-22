#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 09:58:49 2021

@author: cghiaus

Import functions for EPW data files.
Adapted from
https://github.com/pvlib/pvlib-python/blob/master/pvlib/iotools/epw.py
"""

import numpy as np
import pandas as pd
import sys
from scipy.linalg import block_diag


def TCAss(TCd, AssX):
    """
    Parameters
    ----------
    TCd : dictionary of thermal circuits
        DESCRIPTION.
        Dictionary of disassembled thermal circuitss. Example:
            TCd = {'0': TCd0,
                   '1': TCd1,
                   ...
                   'n': TCdn}
        Each thermal circuit is a dictionary:
            TCdk = {'A': A, 'G': G, 'b': b, 'C': C, 'f': f, 'y': y}

    AssX : np.array
        DESCRIPTION.
        Assembling matrix:
            [[TC<-, node, <-TC, node],
             ...
             [TC<-, node, <-TC, node]]

    Returns
    -------
    TCa : Dictionary
        DESCRIPTION.
        Assembled thermal circuit:
            TCdk = {'A': A, 'G': G, 'b': b, 'C': C, 'f': f, 'y': y}
    """
    # Create assembing matrix Ass from AssX
    TCdf = pd.DataFrame(TCd).transpose()

    # Global indexes of the 1st node of each TC
    size_f_eachTCd = TCdf.f.apply(lambda x: np.size(x))
    TCdf['global 1st node'] = np.cumsum(size_f_eachTCd)
    TCdf['global 1st node'] = TCdf['global 1st node'].shift(1).fillna(0)

    # Global indexes of the 1st branch of each TC
    size_b_eachTCd = TCdf.b.apply(lambda x: np.size(x))
    TCdf['global 1st branch'] = np.cumsum(size_b_eachTCd)
    TCdf['global 1st branch'] = TCdf['global 1st branch'].shift(1).fillna(0)

    # Ass = global_1st_node of TC from AssX + its local node
    Ass = np.array([TCdf['global 1st node'][AssX[:, 0]] + AssX[:, 1],
                    TCdf['global 1st node'][AssX[:, 2]] + AssX[:, 3]])
    Ass = Ass.astype(int)

    # Disassembling matrix for temperatures Adθ
    # - matrix that keeps the indexes of themperature nodes
    Adθ = np.eye(sum(size_f_eachTCd))
    # - add the columns that merge
    Adθ[:, Ass[0]] = Adθ[:, Ass[0]] + Adθ[:, Ass[1]]
    # - eliminate the columns that correspond to eliminated nodes
    Adθ = np.delete(Adθ, Ass[1], 1)

    Adq = np.eye(sum(size_b_eachTCd))

    # Ad for [q1 q2 ... θ1 θ2 ...]
    Ad = block_diag(Adq, Adθ)

    # List of indexes for Adq
    row_Adq_local = TCdf.b.apply(lambda x: np.arange(np.size(x)))
    row_Adq_global = row_Adq_local + TCdf['global 1st branch']
    row_Adq_global = [list(x) for x in row_Adq_global]

    # List of indexes for Adθ
    row_Adθ_local = TCdf.f.apply(lambda x: np.arange(np.size(x)))
    row_Adθ_global = row_Adθ_local + TCdf['global 1st node']
    row_Adθ_global += row_Adq_global[-1][-1] + 1
    row_Adθ_global = [list(x) for x in row_Adθ_global]

    row_Ad = list(zip(row_Adq_global, row_Adθ_global))
    row_Ad = [item for sublist in row_Ad for item in sublist]
    row_Ad = [item for sublist in row_Ad for item in sublist]
    row_Ad = [int(i) for i in row_Ad]

    # Ad for [q1 θ1 q2 θ2 ...]
    Ad = Ad[row_Ad, :]

    TCdf['invG'] = TCdf.G.apply(lambda x: np.linalg.inv(x))

    # Block matrices Ki for each TC
    T = TCdf[['A', 'invG', 'C', 'b', 'f', 'y']].copy()
    T['K'] = ""
    Kd = []
    ubf = uby = None
    for k in range(T.shape[0]):
        T['K'][k] = np.block([[T['invG'][k], T['A'][k]],
                              [-T['A'][k].T, T['C'][k]]])
        Kd = block_diag(Kd, T['K'][k])
        ubf = np.block([ubf, T['b'][k], T['f'][k]])
        uby = np.block([uby, T['b'][k], T['y'][k]])
    Kd = np.delete(Kd, obj=0, axis=0)
    ubf = np.delete(ubf, obj=0, axis=0)
    uby = np.delete(uby, obj=0, axis=0)

    Ka = Ad.T @ Kd @ Ad

    # Elements of the assembled circuit
    # total number of flows
    nq = sum(size_b_eachTCd)
    Ga = np.linalg.inv(Ka[:nq, :nq])
    Aa = Ka[:nq, nq:]
    Ca = Ka[nq:, nq:]

    u = Ad.T @ ubf
    ba = u[:nq]
    fa = u[nq:]     # elements of f for merged nodes > 1
    fa[fa.nonzero()] = 1

    u = Ad.T @ uby
    ya = u[nq:]     # elements of f for merged nodes > 1
    ya[ya.nonzero()] = 1

    TCa = {'A': Aa, 'G': Ga, 'b': ba, 'C': Ca, 'f': fa, 'y': ya}

    TCdf['q local'] = row_Adq_local
    TCdf['q global'] = row_Adq_global
    TCdf['θ local'] = row_Adθ_local
    TCdf['θ glob diss'] = row_Adθ_global
    TCdf['θ glob diss'] = TCdf['θ glob diss'].apply(lambda x: np.array(x) - nq)
    return TCa


def tc2ss(A, G, b, C, f, y):
    """
        Parameters
        ----------
        A : TYPE np.array
            adjancecy (TC connection ) matrix:
            #rows = #heat flow rates; #cols = #temperature nodes

        G : TYPE np.array
            square diagonal matrix of conductances
            #rows = #heat flow rates (or resistances)

        b : TYPE np.array
            vector indicating the presence of temperature sources on branches:
                1 for branches with temperature sources, otherwise 0
        C : TYPE np.array
            square diagonal matrix of capacities
        f : TYPE np.array
            vector indicating the presence of flow sources in nodes:
                1 for nodes with heat sources, otherwise 0
        y : TYPE np.array
            vector indicating the temperatures in the outputs:
                1 for output nodes, otherwise 0

        Returns
        -------
        As state matrix in state equation
        Bs input matrix in state equation
        Cs output matrix in observation equation
        Ds input matrix in observation equation
        Idx{1} nodes with capacities
            {2} branches with temp. sources
            {3} nodes with flow sources
            {4} nodes output temperatures

    """

    rC = np.nonzero(np.diag(C))[0]          # rows of non-zero elements in C
    r0 = np.nonzero(np.diag(C) == 0)[0]     # rows of zero elements in C
    # idx_nonzero = {'C': rC,
    #                'b': np.nonzero(b)[0],
    #                'f': np.nonzero(f)[0],
    #                'y': np.nonzero(y)[0]}

    if rC.size == 0:
        sys.exit('Error in dm4bem.tc2ss: capacity C matrix is zero')

    CC = np.diag(C[np.nonzero(C)])
    K = -A.T @ G @ A

    K11 = K[r0, :][:, r0]
    K12 = K[r0, :][:, rC]
    K21 = K[rC, :][:, r0]
    K22 = K[rC, :][:, rC]

    Kb = A.T @ G
    Kb1 = Kb[r0, :]
    Kb2 = Kb[rC, :]

    # State equation
    As = np.linalg.inv(CC) @ (
        -K21 @ np.linalg.inv(K11) @ K12 + K22)
    Bs = np.linalg.inv(CC) @ np.hstack([
        -K21 @ np.linalg.inv(K11) @ Kb1 + Kb2,
        -K21 @ np.linalg.inv(K11),
        np.eye(CC.shape[0])])
    # re-arragne B s in order of f-sources
    # index B for sources [b f0 fC]
    idx_new = np.hstack([np.arange(b.size), b.size + r0, b.size + rC])
    Bs[:, idx_new] = np.array(Bs)
    # indexes of effective inputs [b f]
    inp = np.hstack([np.nonzero(b)[0], A.shape[0] + np.nonzero(f)[0]])
    # extract actual inputs (inputs <> 0)
    Bs = Bs[:, inp]

    # Ds if outputs are all states
    Ds = np.zeros([y[rC].size, np.hstack([b, f]).size])

    # observation equation for outputs that are not states
    Cso = -np.linalg.inv(K11) @ K12
    Dso = -np.linalg.inv(K11) @ np.hstack(
        [Kb1, np.eye(r0.size), np.zeros([r0.size, CC.shape[0]])])

    # observation equation for any output
    Cx = np.zeros([y.size, As.shape[0]])
    Cs = np.diag(y[rC])
    Cx[rC, :] = Cs
    Cx[r0, :] = Cso
    Cs = Cx[np.nonzero(y)[0], :]

    Dx = np.zeros([y.size, np.hstack([b, f]).shape[0]])
    Dx[r0, :] = Dso     # feed-through if no capacity
    Dx[:, idx_new] = np.array(Dx)   # rearange in order of f-sources
    Ds = Dx[np.nonzero(y)[0], :][:, inp]

    return As, Bs, Cs, Ds


# ===========================================================================
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


def read_epw(filename, coerce_year=None):
    '''
    Read an EPW file into a pandas dataframe.

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

    .. [1] `EnergyPlus documentation, Auxiliary Programs
       <https://energyplus.net/documentation>`_
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
