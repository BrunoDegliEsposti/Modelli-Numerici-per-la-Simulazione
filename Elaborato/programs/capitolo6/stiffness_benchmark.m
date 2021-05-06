mu = linspace(1,1000,21);
perfratio = ones(size(mu));

for i = 1:length(mu)
    odefun = @(t,y) [y(2); -y(1)+mu(i)*y(2)*(1-y(1)^2)];
    t0 = 0;
    T = 2*mu(i);
    y0 = [2; 0];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    tic;
    sol1 = ode23(odefun,[t0,T],y0,opts);
    perf1 = toc();
    tic;
    sol2 = ode23s(odefun,[t0,T],y0,opts);
    perf2 = toc();
    perfratio(i) = perf1/perf2;
end

tab = table(mu',perfratio','VariableNames',{'mu','perfratio'});
writetable(tab,'../../tables/capitolo6/van-der-pol-perfratio.dat');
