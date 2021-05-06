f = @(x) x*(x-1)-1;
x0 = 3;
x1 = 2;
tolla = 1e-12;
tollr = 1e-12;
tollf = 0;
nmax = 100;
[x1,fx1,n,errcode,T1] = secanti(f,x0,x1,tolla,tollr,tollf,nmax,true);
writetable(T1,'../../tables/capitolo2/secanti1.dat');
phi = (1+sqrt(5))/2;
M1 = (2*phi-1)^(1-phi);

f = @(x) 1/x - 1e-6;
x0 = 1;
x1 = 2;
tolla = 1e-12;
tollr = 1e-12;
tollf = 0;
nmax = 100;
[x1,fx1,n,errcode,T2] = secanti(f,x0,x1,tolla,tollr,tollf,nmax,true);
writetable(T2,'../../tables/capitolo2/secanti2.dat');

f = @(x) x/(x^2+1);
x0 = 1;
x1 = 2;
tolla = 1e-12;
tollr = 1e-12;
tollf = 0;
nmax = 20;
[x1,fx1,n,errcode,T3] = secanti(f,x0,x1,tolla,tollr,tollf,nmax,true);
writetable(T3,'../../tables/capitolo2/secanti3.dat');

f = @(x) x*x*x + sign(x)*abs(x)^(1/5);
x0 = 10;
x1 = 9;
tolla = 1e-8;
tollr = 1e-8;
tollf = 0;
nmax = 28;
[x1,fx1,n,errcode,T4] = secanti(f,x0,x1,tolla,tollr,tollf,nmax,true);
writetable(T4,'../../tables/capitolo2/secanti4.dat');

f = @(x) x^2;
x0 = 1;
x1 = 0.9;
tolla = 1e-6;
tollr = 0;
tollf = 0;
nmax = 100;
[x1,fx1,n,errcode,T5] = secanti(f,x0,x1,tolla,tollr,tollf,nmax,true);
writetable(T5,'../../tables/capitolo2/secanti5.dat');

f = @(x) expm1(x)-0.5*x^2;
x0 = 1;
x1 = 0.9;
tolla = 1e-12;
tollr = 1e-12;
tollf = 0;
nmax = 100;
[x1,fx1,n,errcode,T6] = secanti(f,x0,x1,tolla,tollr,tollf,nmax,true);
writetable(T6,'../../tables/capitolo2/secanti6.dat');



