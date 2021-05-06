addpath('../capitolo3');
t0 = 0;
T = 20*pi;
h = logspace(-1,-3,10)';
epsil = 0.6;
y0 = [1-epsil, 0, 0, sqrt((1+epsil)/(1-epsil))];

for s = 1:2
    for k = s:4
        err = zeros(size(h));
        for i = 1:length(h)
            N = round((T-t0)/h(i));
            [y,t] = hbvm_nocorH(0, @kepler, y0, k, s, h(i), N, 0, 1);
            H = kepler(y);
            err(i) = max(abs(H+0.5));
        end
        tab = table(h,err,'VariableNames',{'h','errH'});
        tab_filename = sprintf('../../tables/capitolo8/errH-HBVM%d%d.dat',k,s);
        writetable(tab,tab_filename);
    end
end
