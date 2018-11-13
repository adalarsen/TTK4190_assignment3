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
tstop=15000;        % Sim stop time
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)
                
p0=[1000 700]';      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=60*pi/180;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=0;                % Current on (1)/off (0)
%dc= 5*pi/180;

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
kd_heading = 0; %2*zeta*omegan*(m+km)-d
ki_heading = omegan/10*kp_heading; 


kp_surge = 100;%(m+km)*omegan^2-k;
kd_surge = 0; %2*zeta*omegan*(m+km)-d
ki_surge = 0;%omegan/10*kp; 


psi_d.time = tstart:tsamp:tstop';
psi_d.signals.values = 60+0*sin(0.004*psi_d.time)';
r_d = 0.4*cos(0.004*psi_d.time)*0.004;
r_d = r_d';

load('WP.mat')

sim MSFartoystyring_2 % The measurements from the simulink model are automatically written to the workspace.


%% Path plotting


dec = 10;   % Factor for reducing the amount of data in the plot
track = 0;  % 1 for task 2.7 otherwise 0
%pathplotter(x, y,  psi, tsamp, dec, tstart, tstop, track, WP)
pathplotter(p(:,1),p(:,2),psi,tsamp,dec,tstart,tstop,track,WP);





