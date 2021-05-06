% Dataset https://sparse.tamu.edu/LAW/cnr-2000
% Etichette http://law.di.unimi.it/webdata/cnr-2000/
load('cnr-2000.mat');
L = Problem.A';
[m,~] = size(L);
labels = readlines('cnr-2000.urls');
e = ones(m,1);
v = e/m;
f = sum(L)';
d = (f==0);
F = zeros(m,1);
F(f~=0) = 1./f(f~=0);

% Metodo iterativo per il calcolo del PageRank
p = 0.85;
toll = 1e-6;
nmax = ceil(log(toll)/log(p));
y = v;
for n=1:nmax
    y = L*(F.*y) + v*dot(d,y);
    y = p*y + (1-p)*v;
end

% Estrazione degli URL più importanti
k = 20;
[~,idx] = maxk(y,k);
T = table((1:k)', labels(idx), idx, y(idx), ...
    'VariableNames', {'rank','url','i','y'});
writetable(T,'../../tables/capitolo5/pagerank.dat');
