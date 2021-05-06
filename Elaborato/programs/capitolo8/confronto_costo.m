addpath('../capitolo3');
t0 = 0;
T = 20*pi;
h = 1e-2;
N = round((T-t0)/h);
h = (T-t0)/N;
epsil = 0.6;
y0 = [1-epsil; 0; 0; sqrt((1+epsil)/(1-epsil))];

J = kron([0,1;-1,0],eye(2));
odefun = @(t,y) J*kepler_grad_hamiltonian(y,1,1,1);

kmax = 16;
for s = 1:4
    t_naive = zeros(kmax-s+1,1);
    t_optimized = zeros(kmax-s+1,1);
    for k = s:kmax
        % naive algorithm
        [b,c,A,~,~,~] = RKform(k,s);
        order = 2*s;
        tic();
        y1 = runge_kutta(odefun, t0, T, h, y0, A, b, c, order);
        t_naive(k-s+1) = toc();
        
        % optimized algorithm
        tic();
        y2 = hbvm(0, @kepler, y0', k, s, h, N, 0, -1);
        t_optimized(k-s+1) = toc();
        fprintf('k: %d, s: %d, err: %e\n',k,s,norm(y2(2,:)'-y1));
    end
    tab = table((s:kmax)',t_naive,'VariableNames',{'k','tcpu'});
    tab_filename = sprintf(...
        '../../tables/capitolo8/benchmark-naive-HBVMk%d.dat',s);
    writetable(tab,tab_filename);
    tab = table((s:kmax)',t_optimized,'VariableNames',{'k','tcpu'});
    tab_filename = sprintf(...
        '../../tables/capitolo8/benchmark-optimized-HBVMk%d.dat',s);
    writetable(tab,tab_filename);
end
