addpath('../capitolo3');

J = kron([0,1;-1,0],eye(2));
odefun = @(t,y) J*kepler_grad_hamiltonian(y);
t0 = 0;
T = 4*pi;
h = 3*1e-2;
epsil = 0.6;
y0 = [1-epsil; 0; 0; sqrt((1+epsil)/(1-epsil))];

[~,tab] = adams_bashforth(odefun,t0,T,h,y0,2);
tab.H = kepler_hamiltonian(tab{:,2:5}')';
writetable(tab,'../../tables/capitolo8/keplero-adams-bashforth.dat');

[~,tab] = adams_moulton(odefun,t0,T,h,y0,1);
tab.H = kepler_hamiltonian(tab{:,2:5}')';
writetable(tab,'../../tables/capitolo8/keplero-adams-moulton.dat');

[~,tab] = BDF(odefun,t0,T,h,y0,2);
tab.H = kepler_hamiltonian(tab{:,2:5}')';
writetable(tab,'../../tables/capitolo8/keplero-BDF.dat');

k = 5;
s = 1;
N = round((T-t0)/h);
[y,t] = hbvm(0, @kepler, y0', k, s, h, N, 0, 1);
tab = table(t,'VariableNames',{'t'});
tab.y1 = y(:,1);
tab.y2 = y(:,2);
tab.H = kepler_hamiltonian(y')';
writetable(tab,'../../tables/capitolo8/keplero-HBVM51.dat');
