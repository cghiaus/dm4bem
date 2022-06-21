#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 11:04:48 2021

@author: cghiaus

J.A. Duffie, W. A. Beckman (2013) Solar Engineering of Thermal Processes

Data = [month day hour minute HorizRad]
month   number of the month, from 1 to 12
day     day in the month, from 1 to 31 
hour    hour in the 24 hour system (e.g. 13 for 1 p.m.)
% minute  minute from 0 to 59

DirNRad direct normal radiation, Wh/m2
% DifHRad  diffuse horizontal radiation, Wh/m2
%
% B       slope (tilt) angle in deg: [0 180]; 90-vertical; >90-downward facing
% L       local latitude in deg: [-90 90], north positive
% Z       surface azimuth in deg: [-180 180]; 0-south; west-positive
% albedo  albedo of ground = 0.2 Th-CE 2005 pg. 50
%
% Outputs
% PhiDir  direct radiation on the surface in Wh/m2
% PhiDif  diffuse radiation on the surfaca in Wh/m2
% PhiRef  reflected rariation in Wh/m2
"""


def sol_rad_tilt_surf(weather_data, surface_orientation, albedo):
    