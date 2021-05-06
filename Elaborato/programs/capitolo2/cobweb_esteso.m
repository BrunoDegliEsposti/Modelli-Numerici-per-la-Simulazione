d0 = 10; s0 = 0; a = 4; b = 1; rho = 1; nmax = 20;
p = zeros(nmax,1);
p(1) = 3;
p(2) = (d0-s0-b*p(1))/a;
for n = 3:nmax
    p(n) = d0-s0 - b*(1+rho)*p(n-1) + b*rho*p(n-2);
    p(n) = p(n)/a;
end
