% definizione dei parametri e caricamento dei dati
tau = 5;
L = 120;
m = L/tau;
T = readtable('popolazione_italiana_2019.xls');
alpha = T.alpha;
beta = T.beta(1:m-1);
A = diag(beta,-1);
A(1,:) = alpha;

% sostituzioni in avanti (2019 -> 2099)
nmax = 16;
Y = zeros(m,nmax+1);
Y(:,1) = T.y;
for n = 1:nmax
    Y(:,n+1) = A*Y(:,n);
end

% confronto con autovettore dominante
[v0,lambda0] = eigs(A,1);
v0 = v0 ./ sum(v0);
Ynormalized = Y ./ sum(Y);
