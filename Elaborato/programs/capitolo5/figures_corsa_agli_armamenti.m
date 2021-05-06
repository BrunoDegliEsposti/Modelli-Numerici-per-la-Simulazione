% Caso stabile
A = [-1.2,1; 1,-1.2];
xi_eta = [1; 2];
odefun = @(t,y) A*y + xi_eta;
t0 = 0;
T = 15;
h = 0.1;
y0 = [0.5;5];
[~,tab] = trapezoidal_rule(odefun,t0,T,h,y0);
writetable(tab,'../../tables/capitolo5/corsa-agli-armamenti-1.dat');

% Caso instabile
A = [-0.9,2; 1,-2];
xi_eta = [1; 2];
odefun = @(t,y) A*y + xi_eta;
t0 = 0;
T = 10;
h = 0.1;
y0 = [0.5;5];
[~,tab] = trapezoidal_rule(odefun,t0,T,h,y0);
writetable(tab,'../../tables/capitolo5/corsa-agli-armamenti-2.dat');
