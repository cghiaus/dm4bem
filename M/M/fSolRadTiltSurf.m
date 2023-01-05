function [PhiDir, PhiDif, PhiRef]...
 = fSolRadTiltSurf(month, day, hour, minute, RadNDir, RadHDif, B, Z, L, albedo)
% J.A. Duffie, W. A. Beckman (2013) Solar Engineering of Thermal Processes
%
% Data = [month day hour minute HorizRad]
% month   number of the month, from 1 to 12
% day     day in the month, from 1 to 31 
% hour    hour in the 24 hour system (e.g. 13 for 1 p.m.)
% minute  minute from 0 to 59
% DirNRad direct normal radiation, Wh/m2
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

B = B*pi/180; % slope
Z = Z*pi/180; % azimuth
L = L*pi/180; % latitude

n = datenum(0, month, day); % day number in the year

declination_angle=23.45*sin(360*(284+n)/365*pi/180); % eq. 1.6.1a 
d=declination_angle*pi/180;

hour_angle=((hour+minute/60)-12)*15; % Example 1.6.1
h=hour_angle*pi/180;

theta = acos(sin(d)*sin(L)*cos(B) - sin(d)*cos(L)*sin(B)*cos(Z) ...
  + cos(d)*cos(L)*cos(B).*cos(h) + cos(d)*sin(L)*sin(B)*cos(Z).*cos(h)...
  + cos(d)*sin(B)*sin(Z).*sin(h)); % incidence angle eq. 1.6.2
  
theta(theta>pi/2) = pi/2;
PhiDir = RadNDir.*cos(theta); %Th-CE 2005 Eq. 
PhiDir(PhiDir<0) = 0;

PhiDif = RadHDif.*(1 + cos(theta))/2; %Th-CE 2005 Eq. 79

gamma = asin(cos(d)*cos(L).*cos(h) + sin(d)*sin(L));
gamma(gamma<10e-5) = 10e-5;
RadDh = RadNDir.*sin(gamma);
% radiation reflected by the ground:
PhiRef = (RadDh + RadHDif)*albedo*(1 - cos(B)/2);
% endfunction