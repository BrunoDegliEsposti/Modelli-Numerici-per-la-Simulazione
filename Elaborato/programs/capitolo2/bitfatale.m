nmax = 81;
y = zeros(nmax,1);
y(1) = 1;
y(2) = (1-sqrt(5))/2;
for n = 3:nmax
    y(n) = y(n-1) + y(n-2);
end
T = table((0:nmax-1)', abs(y),'VariableNames',{'n','absyn'});
writetable(T,'../../tables/capitolo2/bitfatale.dat');
