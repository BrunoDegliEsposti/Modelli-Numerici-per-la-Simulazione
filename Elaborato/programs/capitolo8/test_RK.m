odefun = @(t,y) [0,-1; 1,0]*y;
t0 = 0;
T = 2*pi;
y0 = [1;0];

[A,b,c,order] = GL2();


n = 7;
e = zeros(n,1);
for i = 1:n
    h = 0.1/(2^(i-1));
    yT = runge_kutta(odefun,t0,T,h,y0,A,b,c,order);
    e(i) = norm(y0-yT);
end
r = e(1:n-1)./e(2:n);
disp(r);
