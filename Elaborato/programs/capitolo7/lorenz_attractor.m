fig1 = figure(1); fig1.Units = 'centimeters';
fig1.InnerPosition(3:4) = [30 30];

sigma = 10;
beta = 8/3;
rho = 28;
odefun = @(t,y) [sigma*(y(2)-y(1)); ...
                 rho*y(1)-y(2)-y(1)*y(3); ...
                 y(1)*y(2)-beta*y(3)];
for i = 1:20
    P = rand(3,1);
    t0 = 0;
    T = 200;
    opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
    sol = ode45(odefun,[t0,T],P,opts);
    t = linspace(T/2,T,10000);
    xyz = deval(sol,t);
    x = xyz(1,:);
    y = xyz(2,:);
    z = xyz(3,:);
    plot3(x,y,z);
    hold on;
end

xlabel('x'); ylabel('y'); zlabel('z'); grid on;
C1 = [+sqrt(beta*(rho-1)), +sqrt(beta*(rho-1)), rho-1];
C2 = [-sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), rho-1];
text(C1(1),C1(2),C1(3),'C1'); text(C2(1),C2(2),C2(3),'C2');
exportgraphics(fig1,'../../figures/capitolo7/attrattore-di-lorenz.png',...
    'ContentType','image','Resolution',300,'BackgroundColor','white');
