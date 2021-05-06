addpath('../capitolo3');

%% Esempio 1
odefun = @(t,y) [-2,1; 1,-2]*y + [1; 2];
t0 = 0;
T = 15;
h = 0.05;
y0 = [0.5;5];
[~,tab] = trapezoidal_rule(odefun,t0,T,h,y0);
plot(tab.t,tab.y1);
hold on;
plot(tab.t,tab.y2);

%% Esempio 2
odefun = @(t,y) [-0.9,2; 1,-2]*y + [1; 2];
t0 = 0;
T = 15;
h = 0.05;
y0 = [0.5;5];
[~,tab] = trapezoidal_rule(odefun,t0,T,h,y0);
plot(tab.t,tab.y1);
hold on;
plot(tab.t,tab.y2);

%% Esempio 3
A = [-1.1,1; 1,-1.1];
xi_eta = [1; 1];
odefun = @(t,y) A*y + xi_eta;
t0 = 0;
T = 55;
h = 0.1;
y0 = [0.5;5];
[~,tab] = trapezoidal_rule(odefun,t0,T,h,y0);
plot(tab.t,tab.y1);
hold on;
plot(tab.t,tab.y2);

%% Esempio 4
A = [-5,1; 1,-5];
xi_eta = [1; 1];
odefun = @(t,y) A*y + xi_eta;
t0 = 0;
T = 3;
h = 0.01;
y0 = [0.5;5];
[~,tab] = trapezoidal_rule(odefun,t0,T,h,y0);
plot(tab.t,tab.y1);
hold on;
plot(tab.t,tab.y2);
