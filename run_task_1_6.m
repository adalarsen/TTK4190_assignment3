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
tstop=10000;        % Sim stop time
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)
                
p0=zeros(2,1);      % Initial position (NED)
v0=[0.001 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=0;                % Current on (1)/off (0)
dc = 0;
nc_arr = [3 7 8.5];

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

result = zeros(3);
for i=1:3
    nc = nc_arr(i);
sim MSFartoystyring_1_6 % The measurements from the simulink model are automatically written to the workspace.

% time-series
tdata = tout;
udata = v(1:end-1,1);   

% nonlinear least-squares parametrization: T dr/dt + u = K delta, delta =
% -dc -> n_c = n_c
% x(1) = 1/T and x(2) = K
x0 = [0.1 1]';
if nc == 3
    F = inline('exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*3','x','tdata');
elseif nc == 7
    F = inline('exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*7','x','tdata');
elseif nc == 8.5
    F = inline('exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*8.5','x','tdata');
end
x = lsqcurvefit(F,x0, tdata, udata);

% estimated parameters
T_surge = 1/x(1);
K_surge = x(2);

result(i,:) = [nc T_surge K_surge];

figure(5)
if nc == 3
    plot(tdata,udata,tdata,exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*3,'--'),grid
elseif nc == 7
    plot(tdata,udata,tdata,exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*7,'--'),grid
elseif nc == 8.5
    plot(tdata,udata,tdata,exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*8.5,'--'),grid
end
title('Least-squares fit of ship-model'),xlabel('time (s)')
legend('Ship n_c = 3','Model, n_c = 3', 'Ship n_c = 7','Model, n_c = 7', 'Ship n_c = 8.5','Model, n_c = 8.5')
hold on;
end
result

% figure(1); clf;
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