alpha = 0.6; rho = 2; % rho maggiore di 1/alpha
G = 800;              % miliardi di euro
Ybar = G/(1-alpha);   % 2000 miliardi
nmax = 20;
Y = zeros(nmax,1);
Y(1) = Ybar-100;
Y(2) = Ybar;
for n = 3:nmax
    Y(n) = alpha*Y(n-1) + rho*alpha*(Y(n-1)-Y(n-2)) + G;
end
T = table((0:nmax-1)', Y, 'VariableNames', {'n','Yn'});
writetable(T,'../../tables/capitolo2/pil.dat');
