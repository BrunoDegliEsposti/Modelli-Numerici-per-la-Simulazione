odefun = @(t,y) [0,-1; 1,0]*y;
t0 = 0;
T = 2*pi;
h = logspace(-1,-3,10);
y0 = [1;0];

% Adams-Bashforth
for k = 1:6
    err = zeros(size(h));
    for i = 1:length(h)
        yT = adams_bashforth(odefun,t0,T,h(i),y0,k);
        err(i) = norm(y0-yT);
    end
    tab = table(h',err','VariableNames',{'h','err'});
    s = sprintf('../../tables/capitolo3/convergenza-adams-bashforth-%d.dat',k);
    writetable(tab,s);
end

% Adams-Moulton
for k = 1:6
    err = zeros(size(h));
    for i = 1:length(h)
        yT = adams_moulton(odefun,t0,T,h(i),y0,k);
        err(i) = norm(y0-yT);
    end
    tab = table(h',err','VariableNames',{'h','err'});
    s = sprintf('../../tables/capitolo3/convergenza-adams-moulton-%d.dat',k);
    writetable(tab,s);
end

% BDF
for k = 1:6
    err = zeros(size(h));
    for i = 1:length(h)
        yT = BDF(odefun,t0,T,h(i),y0,k);
        err(i) = norm(y0-yT);
    end
    tab = table(h',err','VariableNames',{'h','err'});
    s = sprintf('../../tables/capitolo3/convergenza-BDF-%d.dat',k);
    writetable(tab,s);
end
