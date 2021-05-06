mu = linspace(1,1000,21);
sigma = zeros(size(mu));
kappa = zeros(size(mu));
gamma = zeros(size(mu));

for i = 1:length(mu)
    odefun = @(t,y) [y(2); -y(1)+mu(i)*y(2)*(1-y(1)^2)];
    t0 = 0;
    T = 2*mu(i);
    y0 = [2; 0];
    [sigma(i),kappa(i),gamma(i)] = stiffness_ratio(odefun,t0,T,y0);
end

%plot(mu,sigma);
tab = table(mu',sigma','VariableNames',{'mu','sigma'});
writetable(tab,'../../tables/capitolo6/van-der-pol-stiffness.dat');