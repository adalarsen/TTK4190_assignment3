%% Plot 2x2 (debugging)

figure(1); clf;
subplot(2,2,1)
plot(t,psi*180/pi,'b')
hold on
plot(psi_d.time,psi_d*180/pi,'r')
hold on
plot(t,(psi-psi_d)*180/pi,'k')
hold on
legend({'$\psi$','$\psi_d$','$\tilde{\psi}$'},'Interpreter','latex')%,'Location','southeast')
title('Heading (yaw)')
xlabel('Time [s]')
ylabel('Angle [deg]')
set(gca,'FontSize',16)

subplot(2,2,2)
plot(p(:,1),p(:,2),'b')
hold on
%legend({'$\chi$','$\chi_{ref}$'},'Interpreter','latex','Location','southeast')
title('Position')
xlabel('East [m]')
ylabel('North [m]')
set(gca,'FontSize',16)

subplot(2,2,3)
plot(t,r*180/pi,'b')
hold on
plot(t,r_d*180/pi,'r')
hold on
plot(t,(r-r_d)*180/pi,'k')
hold on
legend({'$r$','$r_d$','$\tilde{r}$'},'Interpreter','latex')
title('Yaw rate')
ylabel('Angular rate [deg/s]')
xlabel('Time [s]')
set(gca,'FontSize',16)

subplot(2,2,4)
plot(t,v(:,1),'b')
hold on
plot(t,v(:,2),'r')
hold on

legend({'$u$','$v$'},'Interpreter','latex')
title('Velocity')
ylabel('Aasd [deg]')
xlabel('Time [s]')
set(gca,'FontSize',16)

%% Plot 1.4

figure(4); clf;
subplot(3,1,1)
plot(t,psi*180/pi,'b')
hold on
plot(t,psi_d*180/pi,'r')
hold on
plot(t,(psi-psi_d)*180/pi,'k')
hold on
legend({'$\psi$','$\psi_d$','$\tilde{\psi}$'},'Interpreter','latex')%,'Location','southeast')
title('Heading (yaw)')
ylabel('Angle [deg]')
set(gca,'FontSize',16)
ylim([-30 30])

subplot(3,1,2)
plot(t,r*180/pi,'b')
hold on
plot(t,r_d*180/pi,'r')
hold on
plot(t,(r-r_d)*180/pi,'k')
hold on
legend({'$r$','$r_d$','$\tilde{r}$'},'Interpreter','latex')
title('Yaw rate')
ylabel('Angular rate [deg/s]')
set(gca,'FontSize',16)

subplot(3,1,3)
plot(t,dc*180/pi,'b')
hold on
plot(t,ones(1,length(t))*25,'k--')
plot(t,ones(1,length(t))*-25,'k--')
legend({'$\delta_c$','Saturation limits'},'Interpreter','latex')
title('Rudder input')
ylabel('Angle [deg]')
xlabel('Time [s]')
set(gca,'FontSize',16)
ylim([-32 32])

%% Plot 1.8

figure(3); clf;
subplot(3,1,1)
plot(t,v(:,1),'b')
hold on
plot(t,u_r,'b--')
hold on
plot(t,v(:,1)-u_d,'k')
legend({'$u$', '$u_r$', '$\tilde{u}$'},'Interpreter','latex')%,'Location','southeast')
title('Surge speed')
ylabel('Speed [m/s]')
set(gca,'FontSize',16)
ylim([-5 10])

subplot(3,1,2)
yyaxis left
plot(t,psi*180/pi,'b')
ylabel('Angle [deg]')
lim = max(abs(psi*180/pi));
ylim([-lim lim])
hold on
yyaxis right
plot(t,r*180/pi,'r')
lim = max(abs(r*180/pi));
ylim([-lim lim])
legend({'$\psi$','$r$'},'Interpreter','latex')
title('Heading and yaw rate')
ylabel('Angular rate [deg/s]')
set(gca,'FontSize',16)

subplot(3,1,3)
yyaxis left
plot(t,dc*180/pi,'b')
hold on
plot(t,ones(1,length(t))*25,'b--','HandleVisibility','off')
plot(t,ones(1,length(t))*-25,'b--','HandleVisibility','off')
ylim([-26 26])
ylabel('Angle [deg]')
yyaxis right
plot(t,nc_sat,'r')
hold on
plot(t,ones(1,length(t))*nc_max,'r--','HandleVisibility','off')
plot(t,ones(1,length(t))*0,'r--','HandleVisibility','off')
lim = nc_max +1;
ylim([-lim lim])
legend({'$\delta_c$','$n_c$'},'Interpreter','latex')
title('Control input')
ylabel('[rad/s]')
xlabel('Time [s]')
set(gca,'FontSize',16)

%% Plot 2.4

figure(24); clf;
plot((mod(chi+pi,2*pi)-pi)*180/pi,'b')
hold on
plot((mod(psi+pi,2*pi)-pi)*180/pi,'k')
hold on
plot((mod(chi_d+pi,2*pi)-pi)*180/pi,'r')
hold on
plot(beta*180/pi,'g')
plot(-180*ones(620,1),'k--')
plot(180*ones(620,1),'k--')
plot(0*ones(620,1),'k--')

legend({'$\chi$','$\psi$','$\chi_d$','$\beta$'},'Interpreter','latex')%,'Location','southeast')
title('Course, heading, desired course and crab angle')
ylabel('Angle [deg]')
set(gca,'FontSize',16)
ylim([-200 200])
xlim([0 620])
