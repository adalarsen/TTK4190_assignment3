clear all;
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
tstop=5000;        % Sim stop time
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)

dc_max = 25*(pi/180);
nc_max = (85/60)*2*pi;
p0=zeros(2,1);      % Initial position (NED)
v0=[3 0]';          % Initial velocity (body)
surge_step = [v0(1) 7];
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=1;                % Current on (1)/off (0)

%% heading controller parameters
K= -0.0594;
T= 122.6001;

omegab = 0.05;
zeta = 0.8;
omegan = sqrt(1/(1-2*zeta^2 + sqrt(4*zeta^4-4*zeta^2+2)))*omegab;

m= T/K;
d = 1/K;
k = 0;
km = 0; %optional acceleration feedback

kp_heading = (m+km)*omegan^2-k;
kd_heading = 2*zeta*omegan*(m+km)-d;
ki_heading = omegan/10*kp_heading;  

%% surge controller parameters
d_1 = 0;
d_2 = 1.0375;
m_u = 5500;

zeta_surge = 0.8;
omega_surge = 0.01;
lambda = 0.005;

kp_surge = 2*lambda;
ki_surge = lambda^2;


psi_d.time = tstart:tsamp:tstop';
psi_d.signals.values = 0*psi_d.time';

u_r_0 = 3;
u_r.time = tstart:tsamp:tstop';
u_r.signals.values = 0*u_r.time';
u_r.signals.values(1:124) = surge_step(1)*ones(124,1)';
u_r.signals.values(125:(tstop/tsamp + 1)) = surge_step(2)*ones((length(u_r.time)-124),1)';

u_r_plot = u_r.signals.values;

sim MSFartoystyring_1_8_new % The measurements from the simulink model are automatically written to the workspace.
%% Plot

figure(2)
plot(t,v(:,1),'b')
hold on
plot(t,u_r,'r--')
hold on
plot(t,u_d,'r')
hold on
plot(t,v(:,2),'b--')
hold on
plot(t,nc_sat,'g')
hold on
plot(t,u_d_dot,'g')

legend({'$u$', '$u_r$','$u_d$','$v$','$n_c$', '$\dot{u}_d$'},'Interpreter','latex')
title('Velocity')
ylabel('Speed [m/s] / Shaft speed [rad/s]')
xlabel('Time [s]')
set(gca,'FontSize',16)