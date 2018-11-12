%% Plot 2x2 (debugging)

figure(1); clf;
subplot(2,2,1)
plot(t,psi*180/pi,'b')
hold on
plot(psi_d.time,psi_d.signals.values*180/pi,'r')
hold on
plot(t,(psi-psi_d.signals.values)*180/pi,'k')
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

figure(2); clf;
subplot(3,1,1)
plot(t,psi*180/pi,'b')
hold on
plot(psi_d.time,psi_d.signals.values*180/pi,'r')
hold on
plot(t,(psi-psi_d.signals.values)*180/pi,'k')
hold on
legend({'$\psi$','$\psi_d$','$\tilde{\psi}$'},'Interpreter','latex')%,'Location','southeast')
title('Heading (yaw)')
ylabel('Angle [deg]')
set(gca,'FontSize',16)

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
plot(t,delta_c*180/pi,'b')
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
legend({'$u$'},'Interpreter','latex')%,'Location','southeast')
title('Surge speed')
ylabel('Speed [m/s]')
set(gca,'FontSize',16)

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
plot(t,delta_c*180/pi,'b')
hold on
plot(t,ones(1,length(t))*25,'b--')
plot(t,ones(1,length(t))*-25,'b--')
lim = max(abs(delta_c*180/pi));
ylim([-lim lim])
ylabel('Angle [deg]')
yyaxis right
plot(t,n_c*180/pi,'r')
hold on
plot(t,ones(1,length(t))*(85*2*180/60),'r--')
plot(t,ones(1,length(t))*0,'r--')
lim = max(abs(n_c*180/pi));
ylim([-lim lim])
legend({'$\delta_c$','$n_c$'},'Interpreter','latex')
title('Rudder input')
ylabel('Angle [deg]')
xlabel('Time [s]')
set(gca,'FontSize',16)
