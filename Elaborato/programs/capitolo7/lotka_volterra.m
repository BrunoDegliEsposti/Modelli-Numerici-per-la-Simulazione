a = 2; b = 0.04;
c = 1; d = 0.001;
y0 = [800;60];
epsil = @(t) -0.2/(1+exp(-t+20));
odefun = @(t,y) [( a-epsil(t))*y(1) - b*y(1)*y(2);...
                 (-c-epsil(t))*y(2) + d*y(1)*y(2)];
t0 = 0;
T = 40;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol = ode45(odefun,[t0,T],y0,opts);
t  = sol.x;
y1 = sol.y(1,:);
y2 = sol.y(2,:);
ratio = y1./y2;
tab = table(t',y1',y2',ratio','VariableNames',{'t','y1','y2','ratio'});
writetable(tab,'../../tables/capitolo7/lotka-volterra.dat');

% figure(1); plot(t,y1); title('Prede');
% figure(2); plot(t,y2); title('Predatori');
% figure(3); plot(t,ratio); title('Prede/Predatori');
%figure(4); plot(y1,y2); title('Piano delle fasi');