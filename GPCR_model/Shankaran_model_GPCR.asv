% Recreating Shankaran model matlab
clc
clear
close all

% Parameters
kon = 8.4e7;
koff = 0.37;
kfr = 10;
krr = 10;
kds = 0.065;
ka = 1e-7;
ki = 2e-1;
RT = (5.5e4);
GT = (1e5);
V = 4e-10;
Kd = koff/kon;
Nav = 6.022e23;
A = 0.001;

% Inital conditions

R = RT; 
G = GT;
C = 0;
Ca = C;
Cd = C;
Ga = C;
L = 0.01*Kd;
ft = A*Kd*koff;

% ODE components ODE with the local function
time = 0:250:250;

dVdt = @(time,X) [(-kon*X(1)*X(2)) + koff*X(3);
    (-kon*X(1).*X(2) + koff*X(3))*(1/(Nav*V)) + ft*0;
    kon*X(1).*X(2) - koff*X(3) - kfr*X(3) + krr*X(4);
    kfr*X(3) - krr*X(4) - kds*X(4);
    kds*X(4);
    -ka*X(6).*X(4) + ki*X(7);
    ka*X(6).*X(4) - ki*X(7)];

[time,X] = ode15s(dVdt,time,[R;L;C;Ca;Cd;G;Ga]); % calling the local funtion to be used in ODE45

%  Plotting the result of the entire model
figure(1)
plot(time,X,LineWidth=2)
ylim([0 inf])
legend('R','L','C','Ca','Cd','G','Ga')
xlabel('time')
ylabel('Concentation')
title('ODE solver with local funtion call')

% % Plot both on the same graph
% figure(2);subplot(2,1,1)
% xlabel('time'); ylabel('Concentation')
% plot(time,X,LineWidth=2); legend('R','L','C','Ca','Cd','G','Ga')
% subplot(2,1,2)
% plot(time,X(:,7),LineWidth=2); legend('Species of interest')

% Plotting the Ga component of the model
figure(3)
hold on
for i = [1e-2 1e-1 1e-0 1e1 1e2]
    kds = 0.065*i;
    [time,X] = ode15s(dVdt,time,[R;L;C;Ca;Cd;G;Ga]);
    store = (sprintf('kds = %.1f',kds));
    plot(time,X(:,7)/1000,'DisplayName',store,LineWidth=2)
    ylim([0 inf])
end
hold off
legend show
xlabel('time')
ylabel('Concentation')
title('Ga species over time w local function call')


%% Exacuting with user defined funtion

% Inital conditions
L = 0.01*Kd;
ft = A*Kd*koff;

% Plot of Model
[t,Val] = ode45(@(t,Val) dXdt(t,Val,kds),time,[R;L;C;Ca;Cd;G;Ga]);
figure(4)
plot(time,Val,LineWidth=2) 
ylim([0 inf])
legend('R','L','C','Ca','Cd','G','Ga')
xlabel('time')
ylabel('Concentation')

% Using the ODE stiff solver
figure(5)
[t,Val] = ode15s(@(t,Val) dXdt(t,Val,kds) ,time, [R;L;C;Ca;Cd;G;Ga]);
plot(time,Val(:,7)/1000,LineWidth=2)
ylim([0 inf])
legend('Ga')
xlabel('time')
ylabel('Concentation')

%% ODE user defined funtion
function xdot = dXdt(t,X,kds)
%

% parameters
kon = 8.4e7;
koff = 0.37;
kfr = 10;
krr = 10;

ka = 1e-7;
ki = 2e-1;
V = 4e-10;
Kd = koff/kon;
Nav = 6.022e23;

ft = 0;%0.01*Kd;

xdot(1) = (-kon*X(1).*X(2)) + koff*X(3);
xdot(2) = (-kon*X(1).*X(2) + koff*X(3))*(1/(Nav*V)) + ft;
xdot(3) = kon*X(1).*X(2) - koff*X(3) - kfr*X(3) + krr*X(4);
xdot(4) = kfr*X(3) - krr*X(4) - kds*X(4);
xdot(5) = kds*X(4);
xdot(6) = -ka*X(6).*X(4) + ki*X(7);
xdot(7) = ka*X(6).*X(4) - ki*X(7);
xdot = xdot';
end

