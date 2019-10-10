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
tstop=10000;        % Sim stop time
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)
                
p0=zeros(2,1);      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=1;                % Current on (1)/off (0)
nc = 7.3;

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
psi_d.signals.values = 0.4*sin(0.004*psi_d.time)';
r_des = 0.4*cos(0.004*psi_d.time)*0.004;
r_des = r_des';

sim MSFartoystyring_1_4 % The measurements from the simulink model are automatically written to the workspace.
%% Plot

% figure;
% subplot(2,2,1)
% plot(t,psi*180/pi,'b')
% hold on
% plot(psi_d.time,psi_d.signals.values*180/pi,'r')
% hold on
% plot(t,(psi-psi_d.signals.values)*180/pi,'k')
% hold on
% legend({'$\psi$','$\psi_d$','$\tilde{\psi}$'},'Interpreter','latex')%,'Location','southeast')
% title('Heading (yaw)')
% xlabel('Time [s]')
% ylabel('Angle [deg]')
% set(gca,'FontSize',16)
% 
% 
% subplot(2,2,2)
% plot(p(:,2),p(:,1),'b')
% hold on
% %legend({'$\chi$','$\chi_{ref}$'},'Interpreter','latex','Location','southeast')
% title('Position')
% xlabel('East [m]')
% ylabel('North [m]')
% set(gca,'FontSize',16)
% 
% subplot(2,2,3)
% plot(t,r*180/pi,'b')
% hold on
% plot(t,r_d*180/pi,'r')
% hold on
% plot(t,(r-r_d)*180/pi,'k')
% hold on
% legend({'$r$','$r_d$','$\tilde{r}$'},'Interpreter','latex')
% title('Yaw rate')
% ylabel('Angular rate [deg/s]')
% xlabel('Time [s]')
% set(gca,'FontSize',16)
% 
% 
% subplot(2,2,4)
% plot(t,v(:,1),'b')
% hold on
% plot(t,v(:,2),'r')
% hold on
% 
% legend({'$u$','$v$'},'Interpreter','latex')
% title('Velocity')
% ylabel('Aasd [deg]')
% xlabel('Time [s]')
% set(gca,'FontSize',16)