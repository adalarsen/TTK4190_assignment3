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
tstop=7000;        % Sim stop time
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)

dc_max = 25*(pi/180);
nc_max = (85/60)*2*pi;
p0=[1000 700]';      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=60*pi/180; %60            % Inital yaw angle
r0=0;               % Inital yaw rate
c=1;                % Current on (1)/off (0)
%dc= 5*pi/180;

%% heading controller parameters
K= -0.0594;
T= 122.6001;

omegab = 0.05;
zeta = 0.6;%1%0.6;
omegan = sqrt(1/(1-2*zeta^2 + sqrt(4*zeta^4-4*zeta^2+2)))*omegab;

m= T/K;
d = 1/K;
k = 0;
km = 0; %optional acceleration feedback

kp_heading = (m+km)*omegan^2-k;
kd_heading = 2*zeta*omegan*(m+km)-d;
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
kd_surge = 0;%20;
ki_surge = 0;%0.1;

%% Sim
psi_d.time = tstart:tsamp:tstop';
psi_d.signals.values = 60+0*sin(0.004*psi_d.time)';

load('WP.mat')

% Target
U_t = 3; 
chi_t = atan2(WP(2,2)-WP(2,1),WP(1,2)-WP(1,1));

kappa = nc_max - U_t;

delta_e = 3*304.8; % Look-ahead distance
delta_s = 3*304.8;
e_d = 2*304.8;
s_d = 1*304.8;


sim MSFartoystyring_Complete % The measurements from the simulink model are automatically written to the workspace.


%% Path plotting

dec = 20;   % Factor for reducing the amount of data in the plot
track = 0;  % 1 for task 2.7 otherwise 0
%pathplotter(x, y,  psi, tsamp, dec, tstart, tstop, track, WP)
pathplotter(p(:,1),p(:,2),psi,tsamp,dec,tstart,tstop,track,WP);

%figure(1); hold on
%scatter(p_t_out(:,2), p_t_out(:,1))
%scatter(y_los, x_los)

%% 2.4

beta = asin(v(:,2)./sqrt(v(:,1).^2+v(:,2).^2));

chi = beta + psi;

%% 2.7

% figure(27); clf;
% plot(p_t(:,2),p_t(:,1))
% title('Position of KNM Helge Ingstad')
% xlabel('North [m]')
% ylabel('East [m]')
% xlim([-2500 1100])
% ylim([-22000 8000])

