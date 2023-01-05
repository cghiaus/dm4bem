clear all
clc %clear console
% Load weather data
fileName = 'FRA_Lyon.csv';
from = 1*30*24;   % start time: from 30 Jan.
period = 5*24;   % simulation period: in hours
[Time,Temp,RadNDir,RadHDif,WDir,WSpeed,month,day,hour,minute]...
    = fReadWeather(fileName,from,period);

% Solar radiadion on a tilted surface 
B = 90; % slope (tilt) angle in deg: [0 180]; 90-vertical; >90-downward facing
Z = 0;  % surface azimuth in deg: [-180 180]; 0-south; west-positive
L = 45; % local latitude in deg: [-90 90], north positive
albedo = 0.2; % albedo of ground = 0.2
[PhiDir, PhiDif, PhiRef] = fSolRadTiltSurf(month, day, hour, minute,...
    RadNDir, RadHDif, B, Z, L, albedo);
    
plot(Time/(24*3600), PhiDir,'b'), hold on %direct on surface
% plot(Time/(24*3600), RadNDir,'g') % direct on normal to sun
plot(Time/(24*3600), PhiDif,'r')  % diffusif on surface
% plot(Time/(24*3600), RadHDif,'k') % diffusif on horizontal surface
plot(Time/(24*3600), PhiRef,'m')  % reflected on surface
title('Solar radiation')
xlabel('Time [days]'), ylabel('\Phi [W/m^2]')
legend('\Phi_d_i_r','\Phi_N_d_i_r','\Phi_D_i_f', '\Phi_N_D_i_f')