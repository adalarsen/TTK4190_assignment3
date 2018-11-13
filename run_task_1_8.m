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

nc_max = (85/60)*2*pi;
p0=zeros(2,1);      % Initial position (NED)
v0=[3 0]';       % Initial velocity (body)
surge_step = [v0(1) 7];
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=0;                % Current on (1)/off (0)

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
kd_heading = 0;%2*zeta*omegan*(m+km)-d;
ki_heading = omegan/10*kp_heading;  

%% surge controller parameters
K_surge = 1;
T_surge = 560.4;

omegab_surge = 0.25;
zeta_surge = 1.2;
omegan_surge = sqrt(1/(1-2*zeta_surge^2 + sqrt(4*zeta_surge^4-4*zeta_surge^2+2)))*omegab_surge;

m_surge = T_surge/K_surge;
d_surge = 1/K_surge;
k_surge = 0;
km_surge = 0; %optional acceleration feedback

kp_surge = (m_surge+km_surge)*omegan_surge^2-k_surge;
kd_surge = 0; %2*zeta*omegan*(m+km)-d;
ki_surge = omegan_surge/10*kp_surge; 

kp_surge = 100;
kd_surge = 0;
ki_surge = 0;


psi_d.time = tstart:tsamp:tstop';
psi_d.signals.values = 0*psi_d.time';
r_d = 0*psi_d.time;
r_d = r_d';

% psi_d.signals.values = 0.4*sin(0.004*psi_d.time)';
% r_d = 0.4*cos(0.004*psi_d.time)*0.004;
% r_d = r_d';

u_d.time = tstart:tsamp:tstop';
u_d.signals.values = 0*u_d.time';
u_d.signals.values(1:124) = surge_step(1)*ones(124,1)';
u_d.signals.values(125:(tstop/tsamp + 1)) = surge_step(2)*ones((length(u_d.time)-124),1)';




u_d_plot = u_d.signals.values;



sim MSFartoystyring_1_8 % The measurements from the simulink model are automatically written to the workspace.
%% Plot
%u_d_plot = u_d;


figure(1); clf;
subplot(2,2,1)
plot(t,psi*180/pi,'b')
hold on
plot(psi_d.time,psi_d.signals.values*180/pi,'r--')
hold on
plot(t,(psi-psi_d.signals.values)*180/pi,'k')
hold on
legend({'$\psi$','$\psi_d$','$\tilde{\psi}$'},'Interpreter','latex')%,'Location','southeast')
title('Heading (yaw)')
xlabel('Time [s]')
ylabel('Angle [deg]')
set(gca,'FontSize',16)

subplot(2,2,2)
plot(p(:,2),p(:,1),'b')
hold on
%legend({'$\chi$','$\chi_{ref}$'},'Interpreter','latex','Location','southeast')
title('Position')
xlabel('East [m]')
ylabel('North [m]')
set(gca,'FontSize',16)

subplot(2,2,3)
plot(t,r*180/pi,'b')
hold on
plot(t,r_d*180/pi,'r--')
hold on
plot(t,(r-r_d)*180/pi,'k')
hold on
legend({'$r$','$r_d$','$\tilde{r}$'},'Interpreter','latex')
title('Yaw rate')
ylabel('Angular rate [deg/s]')
xlabel('Time [s]')
set(gca,'FontSize',16)

figure(2)
plot(t,v(:,1),'b')
hold on
plot(t,u_d_plot,'r--')
hold on
plot(t,v(:,2),'b--')
hold on
plot(t,nc,'g')
ylim([0 10])


legend({'$u$', '$u_d$','$v$','$n_c$'},'Interpreter','latex')
title('Velocity')
ylabel('Speed [m/s] / Shaft speed [rad/s]')
xlabel('Time [s]')
set(gca,'FontSize',16)