% Cube with 2 walls and feed-back
clc, clear all

%Physical values
%****************************
Kp = 1e5;              %P-controller gain: large for precision
Sc = 5*3*3; Si = Sc; Sg = 3*3; %surface [m2]: concrete, insulation, glass
Va = 3*3*3;             %air volume[m3]
rhoa = 1.2; ca = 1000;  %indoor air density; heat capacity
Vpa = 1*Va/3600;        %infiltration and ventilation air: volume/hour

% c: concrete; i: insulation;  g: glass
lamc = 2;       lami = 0.04;    lamg = 1.2;   %[W/m K]
rhoccc = 2.5e6; rhoici = 0.02e6; rhogcg = 2.0e6; %[J/m3 K]
wc = 0.2;       wi = 0.08;      wg = 0.01;    %[m]
epswLW = 0.9;   %long wave wall emmisivity
epswSW = 0.5;   %short wave wall emmisivity

epsgLW = 0.6;   %long wave glass emmisivity
taugSW = 0.7;   %short wave glass transmitance
alphagSW = 0.1; %short wave glass absortivity

sigma = 5.67e-8;%[W/m2 K4]
Fwg = 1/5;      %view factor wall - glass
Tm = 20 + 273;  %mean temp for radiative exchange

% convection coefficients
ho = 10;    hi = 4;   %[W/m2 K]
%***************************

% Conductances and capacities
Gc = lamc/wc*Sc; Cc = Sc*wc*rhoccc; %concrete
Gi = lami/wi*Si; Ci = Si*wi*rhoici; %insulation
Gg = lamg/wg*Sg; Cg = Sg*wg*rhogcg; %glass
Ca = Va*rhoa*ca;
% Convection
Gwo = ho*Sc; Gwi = hi*Si;           %convection wall out; wall in
Ggo = ho*Sg; Ggi = hi*Sg;           %convection glass out; glass in
% Long wave radiative exchange
GLW1 = epswLW/(1-epswLW)*Si*4*sigma*Tm^3;
GLW2 = Fwg*Si*4*sigma*Tm^3;
GLW3 = epsgLW/(1-epsgLW)*Sg*4*sigma*Tm^3;
GLW = 1/(1/GLW1 + 1/GLW2 +1/GLW3);  %long-wave exg. wall-glass
% ventilation & advection
Gv = Vpa*rhoa*ca;                   %air ventilation
% glass: convection outdoor & conduction
Ggs = 1/(1/Ggo + 1/(2*Gg));         %cv+cd glass

% Thermal network
% *****************************************************************
A(1,1) = 1;
A(2,1) = -1; A(2,2) = 1;
A(3,2) = -1; A(3,3) = 1;
A(4,3) = -1; A(4,4) = 1;
A(5,4) = -1; A(5,5) = 1;
A(6,5) = -1; A(6,6) = 1;
A(7,5) = -1; A(7,7) = 1;
A(8,6) = -1; A(8,7) = 1;
A(9,8) = 1;
A(10,6) = 1; A(10,8) = -1;
A(11,7) = 1;
A(12,7) = 1;

G = diag([Gwo 2*Gc 2*Gc 2*Gi 2*Gi GLW Gwi Ggi Ggs 2*Gg Gv Kp]');
b = zeros(12,1); b(1) = 1; b(9) = 1; b(11) = 1; b(12) = 1;
C = diag([0 Cc 0 Ci 0 0 Ca Cg]);
C = diag([0 Cc 0 Ci 0 0 0   0]);
f = zeros(8,1); f(1) = 1; f(5) = 1; f(7) = 1; f(8) = 1;
y = zeros(8,1); 
y(7) = 1;

% Thermal circuit -> state-space
% *********************************************************************
[As,Bs,Cs,Ds] = fTC2SS(A,G,b,C,f,y);

% Maximum time-step
dtmax = min(-2./eig(As))% [s]
dt = 5                  % [s] time step

% Step response
% *********************************************************************
duration = 3600*24*3;   % [s] time duration 
n = floor(duration/dt); % no of time samples

Time = 0:dt:(n-1)*dt;   % time
nth = size(As,1);       % no of state variables
th = zeros(nth,n);      % zero initial conditions
u = zeros(8,n);         % u = [To To To Tsp Phio Phii Qaux Phia]
u(1:3,:) = ones(3,n);   % To = step variation

for k = 1:n-1
 th(:,k+1) = (eye(nth) + dt*As)*th(:,k) + dt*Bs*u(:,k);
end
y = Cs*th + Ds*u;
subplot(3,1,1)
plot(Time/3600,y) 
xlabel('Time [h]'),ylabel('T [C]')
title('Step response for T_o = 1 C'), 

% Simulation with weather data
% Load weather data
fileName = 'FRA_Lyon.csv';
from = 6*30*24 + 25*24;   % start time: from 24 Jan.
period = 10*24; % simulatio, period: for 10 days

[Time,Temp,RadNDir,RadHDif,WDir,WSpeed,month,day,hour,minute]...
    = fReadWeather(fileName,from,period);

B = 90; Z = 0; L = 45; albedo = 0.2;
[PhiDir, PhiDif, PhiRef] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);

% interpolate weather data for time step dt
Temp = interp1(Time, Temp, [Time(1):dt:Time(end)]');
PhiDir = interp1(Time, PhiDir, [Time(1):dt:Time(end)]');
PhiDif = interp1(Time, PhiDif, [Time(1):dt:Time(end)]');
PhiRef = interp1(Time, PhiRef, [Time(1):dt:Time(end)]');
Time = [Time(1):dt:Time(end)]';

n = size(Time,1);
th = zeros(nth,n);
Qa = zeros(n,1);  %auxiliary sources (electrical, persons, etc.)
TintSP = 20*ones(n,1);

% Inputs
PhiTot = PhiDir + PhiDif + PhiRef;
u = [Temp Temp Temp TintSP ...
  epswSW*Sc*PhiTot taugSW*epswSW*Sg*PhiTot Qa alphagSW*Sg*PhiTot]';
% Memory alocation and initial value
th = zeros(size(As,2),n);

for k = 1:n-1
 th(:,k+1) = (eye(nth) + dt*As)*th(:,k) + dt*Bs*u(:,k);
end
y = Cs*th + Ds*u;
subplot(3,1,2)
plot(Time/3600,y,'b',Time/3600,Temp,'g')
xlabel('Time [h]'),ylabel('T [C]')
title('Simulation for weather')

% Solar radiation
subplot(3,1,3)
plot(Time/(24*3600), PhiDir,'b'), hold on %direct on surface
plot(Time/(24*3600), PhiDif,'r')  % diffusif on surface
plot(Time/(24*3600), PhiRef,'m')  % reflected on surface
xlabel('Time [days]'), ylabel('\Phi [W/m^2]')
legend('\Phi_d_i_r','\Phi_D_i_f', '\Phi_R_e_f_D_i_f')
title('Solar radiation')