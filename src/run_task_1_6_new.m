clear;clc;
close all;
%% Information 
% This file is only an example of how you can start the simulation. The
% sampling time decides how often you store states. The execution  time
% will increase if you reduce the sampling time.

% Please note that the file "pathplotter.m" (only used in the second part
% of the assignment) shows the ship path during the path following and
% target tracking part of the assignment. It can be clever to adjust the sampling
% time when you use that file because it draws a sketch of the ship in the
% North-East plane at each time instant. Having a small sampling time will
% lead to multiple ship drawings on top of each other. 

% You should base all of your simulink models on the MSFartoystyring model
% and extend that as you solve the assignment. For your own sake, it is
% wise to create a new model and run file for each task.

% The msfartoystyring.m file includes the ship model. You are not allowed
% to change anything within that file. You need to include that file in
% every folder where you have a simulink model based on
% "MSFartoystyring.slx". 

% WP.mat is a set of six waypoints that you need to use in the second part of
% the assignment. The north position is given in the first row and the east
% position in the second row. 


%%
tstart=0;           % Sim start time
tstop=7000;        % Sim stop time
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)
                
p0=zeros(2,1);      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=0;                % Current on (1)/off (0)
dc = 0;

K= -0.0594;
T= 122.6001;

omegab = 0.05;
zeta = 0.8;
omegan = sqrt(1/(1-2*zeta^2 + sqrt(4*zeta^4-4*zeta^2+2)))*omegab;

m= T/K;
d = 1/K;
k = 0;
km = 0; %optional acceleration feedback

kp = (m+km)*omegan^2-k;
kd = 2*zeta*omegan*(m+km)-d;
ki = omegan/10*kp; 

psi_d.time = tstart:tsamp:tstop';
psi_d.signals.values = 0*psi_d.time';

nc_1 = 8;
nc_2 = 4;

figure(1)
nc = nc_1;
sim MSFartoystyring_1_6
plot(t,v(:,1),'b')
hold on
u_1 = v(end,1);
nc = nc_2;
sim MSFartoystyring_1_6
u_2 = v(end,1);
plot(t,v(:,1),'r--')
hold off
grid on;
legend({'$u_{ship}$','$u_{model}$',},'Interpreter','latex')
title('Velocity')
ylabel('[m/s]')
xlabel('Time [t]')
set(gca,'FontSize',16)

U = [u_1 u_1*abs(u_1);
     u_2 u_2*abs(u_2)];

N = [nc_1*abs(nc_1) nc_2*abs(nc_2)]';

D = inv(U)*N
d_1 = D(1);
d_2 = D(2);

%% m_u estimation

tstop = 5000;
m_u = 5500;
nc = 7.3;
v0=[0.001 0]'; 
sim MSFartoystyring_1_6
m_u = 5500;
sim surge_model


figure(2)
plot(t,v(:,1),'b')
hold on
plot(t,u,'r--')
hold off
grid on;

legend({'$u_{ship}$','$u_{model}$',},'Interpreter','latex')
title('Velocity')
ylabel('[m/s]')
xlabel('Time [t]')
set(gca,'FontSize',16)



