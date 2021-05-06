f = @(t,y) [1+y(1)+y(2);1e8+y(1)+y(2)];
t0 = 0;
y0 = [1e-10;-1]; % y0 = [1e10;-1];
numerical_jacobian(f,t0,y0)
numjacobian(f,t0,y0)
