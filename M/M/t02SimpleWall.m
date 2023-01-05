% Simple wall with capacities in all temperature nodes
% Inputs: outdoor temperature, indoor convection heat flow rate (from HVAC))
clc, clear all

% Physical properties
% *******************
Sw = 3*3;   %wall surface [m2]
Va = 3*3*3; %air volume[m3]
% 1: concrete ; 2: insulation ; 3: air
lam1 = 1.4;   lam2 = 0.04; %[W/m K]
rho1c1 = 2.02e6; rho2c2 = 0.02e6; rho3c3 = 1.2e3;  %[J/K m3]
w1 = 0.2;   w2 = 0.08;  %[m] wall width
x1 = 0.05;  x2 = 0.04;  %[m] discretizaton slice width
% convection coefficents
ho = 10; hi = 4;   %[W/m2 K]
% Thermal resistances
% concrete
Rc = w1/(lam1*Sw);   Cc = Sw*w1*rho1c1;
% insulation
Ri = w2/(lam2*Sw);   Ci = Sw*w2*rho2c2;
% convection
Rvi = 1/(hi*Sw); Rvo = 1/(ho*Sw);

% Dynamic model
%****************
% Thermal circuit
nth = 7; nq = 7;  % # of temperature node, # of flow nodes
% resistances
R = zeros(nq,nq);
R(1,1) = Rvo + Rc/8; 
R(2,2) = Rc/4; R(3,3)=R(2,2); R(4,4)=R(2,2);
R(5,5) = Rc/8 + Ri/4; 
R(6,6) = Ri/2; 
R(7,7) = Ri/4 + Rvi;
G = inv(R);
% capacitances
C = zeros(nth,nth);
C(1,1) = 1/4*Sw*w1*rho1c1; C(2,2)=C(1,1); C(3,3)=C(1,1); C(4,4)=C(1,1);
C(5,5) = 1/2*Sw*w2*rho2c2; C(6,6)=C(5,5);
C(7,7) = Va*rho3c3;
% arc-node incidence matrix
A = eye(nq+1,nth);
A = -diff(A,1,1)';

b = [1 0 0 0 0 0 0]'; f = [0 0 0 0 0 0 0]';
thSteadyTo = inv(A'*G*A)*(A'*G*b + f) % steady-state temperature: input To = 1
b = [0 0 0 0 0 0 0]'; f = [0 0 0 0 0 0 1]'; % idem, for Qh = 1 
thSteadyQh = inv(A'*G*A)*(A'*G*b + f)

% State-space representation
B = inv(C)*[A'*G eye(nth,nth)]; % inputs u = [b; f] size(b)=nq, size(f)=nth;
B = B(:,[1 14]);        % select the 2 relevant inputs: 1->To and 14->Qh
A = inv(C)*(-A'*G*A); 
C = zeros(1,7);C(7)=1;  % output: th(7)
D = [0 0];

% Time integration using Euler forward
% ************************************
% Step response
disp(['max dt = ',num2str(min(-2./eig(A))),'[s]'])
dt = 360
n = 3600/dt*24*30; % no of time samples for 30 days
n = floor(n);
Time = 0:dt:(n-1)*dt; % time
u = [ones(1,n); zeros(1,n)];
th = zeros(nth,n); thi = zeros(nth,n);
for k = 1:n-1
 th(:,k+1) = (eye(nth) + dt*A)*th(:,k) + dt*B*u(:,k);
 thi(:,k+1) = inv((eye(nth) - dt*A))*(thi(:,k) + dt*B*u(:,k));
end
subplot(2,2,1)
plot(Time/3600,th(7,1:n),'r', Time/3600,thi(7,1:n),'b')
xlabel('Time [hours]'), ylabel('T [C]')
title('Step response for T_o = 1 C')

u = [zeros(1,n); ones(1,n)];
for k = 1:n-1
 th(:,k+1) = (eye(nth) + dt*A)*th(:,k) + dt*B*u(:,k);
 thi(:,k+1) = inv((eye(nth) - dt*A))*(thi(:,k) + dt*B*u(:,k));
end
subplot(2,2,2)
plot(Time/3600,th(7,1:n),'r', Time/3600,thi(7,1:n),'b')
xlabel('Time [h]'), ylabel('T [C]')
title('Step response for Q_h = 1 W')

% Simulation with outoor temperature
% Load weather data
fileName = 'FRA_Lyon.csv';
from = 6*30*24 + 25*24;   % start time
period = 30*24; % simulation period

[Time,Temp,RadNDir,RadHDif,WDir,WSpeed,month,day,hour,minute]...
    = fReadWeather(fileName,from,period);

b = 90; z = 0; l = 45; albedo = 0.2;
[PhiDir, PhiDif, PhiRef] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, b, z, l, albedo);

Temp = interp1(Time, Temp, [Time(1):dt:Time(end)]'); %interpolate for dt
Time = [Time(1):dt:Time(end)]';

n = size(Time,1);
th = zeros(nth,n);thi = zeros(nth,n);
Qh = zeros(n,1);
u = [Temp'; Qh'];
for k = 1:n-1
 th(:,k+1) = (eye(nth) + dt*A)*th(:,k) + dt*B*u(:,k);
 thi(:,k+1) = inv((eye(nth) - dt*A))*(thi(:,k) + dt*B*u(:,k));
end

subplot(2,2,3)
plot(Time/3600/24,th(7,1:n),'r', Time/3600/24, Temp,'b')
xlabel('Time [days]'), ylabel('T [C]')
title('Simulation explicit Euler')

subplot(2,2,4)
plot(Time/3600/24,thi(7,1:n),'r', Time/3600/24, Temp,'b')
xlabel('Time [days]'), ylabel('T [C]')
title('Simulation implicit Euler')