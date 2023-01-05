% Cube with 2 walls: assembly
clc, clear all

% Input data
%****************************
% Parameters
Kp = 1e1;               % P-controller gain: large for precision
nc = 4;                 % number of concrete meshes
ni = 3;                 % number of insulation meshes
dt = 5                  % [s] simulation time step

% Physical values
Sc = 5*3*3; Si = Sc;    % surface : concrete, insulation [m2]
Sg = 3*3;               % surface glass
Va = 3*3*3;             % air volume [m3]
rhoa = 1.2; ca = 1000;  % indoor air density; heat capacity
ACH = 3;                % air changes per hour
Vda = ACH*Va/3600;        % infiltration and ventilation air: volume/hour

% c: concrete; i: insulation;  g: glass
lamc = 2;       lami = 0.04;    lamg = 1.2;     % conductivity [W/m K]
rhoccc = 2.5e6; rhoici = 2.0e6; rhogcg = 2.0e6; % [J/m3 K]
wc = 0.2;       wi = 0.08;      wg = 0.01;      % width [m]
epswLW = 0.8;   % long wave wall emmisivity
epswSW = 0.8;   % short wave wall emmisivity

epsgLW = 0.9;   % long wave glass emmisivity
taugSW = 0.7;   % short wave glass transmitance
alphagSW = 0.2; % short wave glass absortivity

sigma = 5.67e-8;% [W/m2 K4]
Fwg = 1/5;      % view factor wall - glass
Tm = 20 + 273;  % mean temp for radiative exchange


hi = 4; ho = 10;% convection coefficients [W/m2 K]
%***************************

% Conductances and capacities
Gc = lamc/wc*Sc; Cs = Sc*wc*rhoccc; % concrete
Gcm = 2*nc*Gc*ones(1,2*nc);         % meshed concrete
Ccm = Cs/nc*mod(0:2*nc-1,2);

Gi = lami/wi*Si; Ci = Si*wi*rhoici; % insulation
Gim = 2*ni*Gi*ones(1,2*ni);         % meshed insulation
Cim = Ci/ni*mod(0:2*ni-1,2);

Gg = lamg/wg*Sg; Cg = Sg*wg*rhogcg; %glass
Ca = Va*rhoa*ca;
% Ca =0; Cg = 0;  % if air and glass are neglected

% Convection
Gwo = ho*Sc; Gwi = hi*Si;           %convection wall out; wall in
Ggo = ho*Sg; Ggi = hi*Sg;           %convection glass out; glass in
% Long wave radiative exchange
GLW1 = epswLW/(1-epswLW)*Si*4*sigma*Tm^3;
GLW2 = Fwg*Si*4*sigma*Tm^3;
GLW3 = epsgLW/(1-epsgLW)*Sg*4*sigma*Tm^3;
GLW = 1/(1/GLW1 + 1/GLW2 +1/GLW3);  %long-wave exg. wall-glass
% advection
Gv = Vda*rhoa*ca;                   %air ventilation
% glass: convection outdoor & conduction
Ggs = 1/(1/Ggo + 1/(2*Gg));         %cv+cd glass


% Thermal network
% *****************************************************************
% Dissembled circuit
nq = 1+2*(nc+ni); nt = 1+2*(nc+ni); % no of flows & temps
A1 = eye(nq+1,nt);    % (#flows+1, #temp)
A1 = -diff(A1,1,1)';
G1 = diag([Gwo Gcm Gim]');
b1 = zeros(nq,1); b1(1) = 1;
C1 = diag([Ccm Cim 0]);
f1 = zeros(nt,1); f1(1) = 1; f1(end) = 1;
y1 = zeros(nt,1);
TCd{1} = {A1,G1,b1,C1,f1,y1};

A2 = [-1 1 0; -1 0 1; 0 -1 1];
G2 = diag([GLW Gwi Ggi]');
b2 = [0 0 0]';
C2 = diag([0 0 Ca/2]');
f2 = [1 0 1]';
y2 = [0 0 1]';
TCd{2} = {A2,G2,b2,C2,f2,y2};

A3 = [1 0;-1 1];
G3 = diag([Ggs 2*Gg]);
b3 = [1 0]';
C3 = diag([Cg 0]);
f3 = [1 0]';
y3 = [0 0]';
TCd{3} = {A3,G3,b3,C3,f3,y3};

A4(1,1) = 1;
A4(2,1) = 1;
G4 = diag([Gv Kp]);
b4 = [1 1]';
C4 = diag([Ca/2]);
f4 = 1;
y4 = 1;
TCd{4} = {A4,G4,b4,C4,f4,y4};

% Assembling matrix
AssX = [1 nt 2 1;...% TC1#5 <- TC2#1
        2 2 3 2;... % TC2#2 <- TC3#2
        2 3 4 1];   % TC2#3 <- TC4#1

% Call Assembling function
[TCa, Idx] = fTCAssAll(TCd, AssX);
A = TCa{1}; G = TCa{2}; b = TCa{3}; C = TCa{4}; f = TCa{5}; y = TCa{6};

% Thermal model -> state-space
% *********************************************************************
[As,Bs,Cs,Ds] = fTC2SS(A,G,b,C,f,y);    % model with controller

% Maximum time-step
dtmax = min(-2./eig(As))% [s]

% Simulation
% *********************************************************************
% Load weather data
fileName = 'FRA_Lyon.csv';
from = 1*30*24 + 25*24;   % start time: from 24 Jan.
period = 10*24; % simulation period: for 10 days

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
nth = size(As,1);       % no of state variables
Qa = zeros(n,1);  %auxiliary sources (electrical, persons, etc.)
TintSP = 20*ones(n,1);

% Inputs
PhiTot = PhiDir + PhiDif + PhiRef;
u = [Temp Temp Temp TintSP/2 ...
  epswSW*Sc*PhiTot taugSW*epswSW*Sg*PhiTot alphagSW*Sg*PhiTot ...
  Qa]';
% Memory alocation and initial value
% th = 20*ones(size(Ac,2),n);
th = zeros(size(As,2),n);
y = zeros(1,n); y(1) = TintSP(1);
qHVAC = zeros(1,n);

for k = 1:n-1
    if y(k) < TintSP(k)
        u(8,k) = 10000*(TintSP(k) - y(k));
        th(:,k+1) = (eye(nth) + dt*As)*th(:,k) + dt*Bs*u(:,k);
        y(:,k+1) = Cs*th(:,k+1) + Ds*u(:,k+1);
        %qHVAC(:,k+1) = Kp*(TintSP(k+1) - y(k+1));
    else
        u(8,k) = 0;
        th(:,k+1) = (eye(nth) + dt*As)*th(:,k) + dt*Bs*u(:,k);
        y(:,k+1) = Cs*th(:,k+1) + Ds*u(:,k+1);
    end
end
subplot(3,1,1)
plot(Time/3600,y,'g',Time/3600,TintSP,'r', Time/3600,Temp,'b')
xlabel('Time [h]'),ylabel('T [C]')
title('Simulation for weather')

subplot(3,1,2)
plot(Time/3600,u(8,:))
xlabel('Time [h]'),ylabel('Q_h [W]')
title('Q_h')
axis([0 inf -1000 1e3])

% Solar radiation
subplot(3,1,3)
plot(Time/(24*3600), PhiDir,'b'), hold on %direct on surface
plot(Time/(24*3600), PhiDif,'r')  % diffusif on surface
plot(Time/(24*3600), PhiRef,'m')  % reflected on surface
xlabel('Time [days]'), ylabel('\Phi [W/m^2]')
title('Solar radiation')
legend('\Phi_d_i_r','\Phi_D_i_f', '\Phi_R_e_f_D_i_f')