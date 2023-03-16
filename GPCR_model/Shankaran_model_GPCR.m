% Recreating Shankaran model matlab
clc
clear
close all

% Parameters
kon = 8.4*10.^7;
koff = 0.37;
kfr = 10;
krr = 10;
kds = 0.065;
ka = 10.^-6;
ki = 2*10.^-1;
RT = 5.5*10.^4;
GT = 1*10.^5;
V = 4*10.^-10;
Kd = koff/kon;
Nav = 6.022*10. ^23;
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

dVdt = @(time,X) [-kon*X(1)*X(2) + koff*X(3);
    (-kon*X(1)*X(2)/(Nav*V)) + L;
    kon*X(1)*X(2) - koff*X(3) - kfr*X(3) + krr*X(4);
    kfr*X(3) - krr*X(4) - kds*X(4);
    kds*X(3);
    -ka*X(6)*X(4) + ki*X(7);
    ka*X(6)*X(4) - ki*X(7)];

[time,X] = ode45(dVdt,time,[R;0;C;Ca;Cd;G;Ga]); % calling the local funtion to be used in ODE45

%  Plotting the result of the entire model
figure(1)
plot(time,X,LineWidth=2)
legend('R','L','C','Ca','Cd','G','Ga')
xlabel('time')
ylabel('Concentation')

% % Plot both on the same graph
% figure(1);subplot(2,1,1)
% xlabel('time'); ylabel('Concentation')
% plot(time,X,LineWidth=2); legend('R','L','C','Ca','Cd','G','Ga')
% subplot(2,1,2)
% plot(time,X(:,7),LineWidth=2); legend('Species of interest')

% Plotting the component of the model
figure(2)
plot(time,X(:,7),LineWidth=2)
legend('Species of interest')
xlabel('time')
ylabel('Concentation')


%% Exacuting with user defined funtion

% Inital conditions
R = RT;
G = GT;
C = 0;
Ca = C;
Cd = C;
Ga = C;
L = 0.01*Kd;
ft = A*Kd*koff;

[t,Val] = ode45(@dXdt,time,[R;0;C;Ca;Cd;G;Ga]);
figure(3)
plot(time,Val,LineWidth=2)
legend('R','L','C','Ca','Cd','G','Ga')
xlabel('time')
ylabel('Concentation')

%% ODE user defined funtion
function xdot = dXdt(t,X)
% parameters
kon = 8.4*10.^7;
koff = 0.37;
Kd = koff/kon;
kfr = 10;
krr = 10;
kds = 0.065;
ka = 10.^-6;
ki = 2*10.^-1;
V = 4*10.^-10;
Nav = 6.022*10. ^23;

ft = 4.40476190476191e-11;

xdot(1) = (-kon*X(1)*X(2)) + koff*X(3);
xdot(2) = (-kon*X(1).*X(2)/(Nav*V)) + ft;
xdot(3) = kon*X(1).*X(2) - koff*X(3) - kfr*X(3) + krr*X(4);
xdot(4) = kfr*X(3) - krr*X(4) - kds*X(4);
xdot(5) = kds*X(3);
xdot(6) = -ka*X(6).*X(4) + ki*X(7);
xdot(7) = ka*X(6).*X(4) - ki*X(7);
xdot = xdot';
end

