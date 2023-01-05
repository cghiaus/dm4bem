function [Time,Temp,RadNDir,RadHDif,WDir,WSpeed,month,day,hour,minute]...
    = fReadWeather(fileName,from,period)
% Reads hourly weather data from .csv file
% Obtain .epw weather data file from https://energyplus.net/weather
% Delete the header (first 8 lines) and the comumn 7 (column F)
% Save the file as .csv
%
% Inputs
% fileName: string, e.g.'FRA_Lyon.csv'
% from: line number in data file from which data are read, e.g. 5*24
% period: incremental vector of the lines which are read, e.g.10*24 
%
% Outputs
% Time: hours 

M = dlmread(fileName);

period = [1:period]';
Time = 3600*period;
Temp = M(from+period,7-1);      %Dry bulb temperarure [°C]
RadNDir =  M(from+period,15-1); %Direct normal solar radiation [Wh/m2]
RadHDif =  M(from+period,16-1); %Diffuse horizontal solar radiation [Wh/m2]
WDir = M(from+period,21-1);     %Wind direcrion: N=00; E=90Â°; S=180Â°; W=270°
WSpeed = M(from+period,22-1);   %wind speed [m/s]

month = M(from+period,2); %
day = M(from+period,3); %
hour = M(from+period,4); %
minute = M(from+period,5); %


% B = 90; Z = 0; L = 45; albedo = 0.2;
% [PhiDir, PhiDif, PhiRef] = fSolRadTiltSurf(month, day, hour, minute, ...
%   RadNDir, RadHDif, B, Z, L, albedo);
% plot(Time/(24*3600), PhiDir,'b'), hold on %direct on surface
% plot(Time/(24*3600), RadNDir,'g') % direct on normal to sun
% plot(Time/(24*3600), PhiDif,'r')  % diffusif on surface
% plot(Time/(24*3600), RadHDif,'k') % diffusif on horizontal surface
% plot(Time/(24*3600), PhiRef,'m')  % reflected on surface
% legend('\Phi_d_i_r','\Phi_N_d_i_r','\Phi_D_i_f', '\Phi_N_D_i_f')
