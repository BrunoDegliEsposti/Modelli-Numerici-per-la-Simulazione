sigma = 10;
beta = 8/3;

%% Origine asintoticamente stabile
rho = 0.5;
for i = 1:100
    y0 = 100*rand(3,1);
    [~,y] = lorenz_solver(0, 100, y0, sigma, beta, rho);
    if norm(y(end,:)) >= 1e-11
        disp(y0);
    end
end

%% Due punti asintoticamente stabili
rho = 1.5;
c1 = [+sqrt(beta*(rho-1)), +sqrt(beta*(rho-1)), rho-1]; % rosso
c2 = [-sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), rho-1]; % blu
plot3(c1(1),c1(2),c1(3),'or',c2(1),c2(2),c2(3),'ob');
hold on;
for i = 1:250
    y0 = 100*(rand(3,1)-0.5);
    [~,y] = lorenz_solver(0, 100, y0, sigma, beta, rho);
    if norm(y(end,:)-c1) <= 1e-8
        plot3(y0(1),y0(2),y0(3),'>r');
    elseif norm(y(end,:)-c2) <= 1e-8
        plot3(y0(1),y0(2),y0(3),'<b');
    else
        plot3(y0(1),y0(2),y0(3),'dg');
    end
end

%% Strano attrattore
rho = 28;
c1 = [+sqrt(beta*(rho-1)), +sqrt(beta*(rho-1)), rho-1]; % rosso
c2 = [-sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), rho-1]; % blu
[t,y] = lorenz_solver(0, 125, [1,1,1], sigma, beta, rho);
plot3(c1(1),c1(2),c1(3),'or',c2(1),c2(2),c2(3),'ob');
hold on;
plot3(y(:,1),y(:,2),y(:,3),'k');







