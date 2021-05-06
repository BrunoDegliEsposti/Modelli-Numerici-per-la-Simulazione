sigma = 10;
beta = 8/3;
rho = 28;
odefun = @(t,y) [sigma*(y(2)-y(1)); ...
                 rho*y(1)-y(2)-y(1)*y(3); ...
                 y(1)*y(2)-beta*y(3)];
P1 = [1;1;1];
P2 = [1+1e-10;1;1];
t0 = 0;
T = 60;
opts = odeset('RelTol',3e-14,'AbsTol',3e-14);
sol1 = ode45(odefun,[t0,T],P1,opts);
sol2 = ode45(odefun,[t0,T],P2,opts);
t = linspace(t0,T,1000);
xyz1 = deval(sol1,t);
xyz2 = deval(sol2,t);
err = vecnorm(xyz1-xyz2);
tab = table(t',err','VariableNames',{'t','err'});
writetable(tab,'../../tables/capitolo7/lorenz.dat');
%semilogy(t,err);
