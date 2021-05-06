odefun = @(t,y) 10*y-2*y*y;
t0 = 0;
T = 2.5;
y0 = 1;
h = 0.01;
[yT,tab] = midpoint_rule(odefun,t0,T,h,y0);
plot(tab.t,tab.y1);
