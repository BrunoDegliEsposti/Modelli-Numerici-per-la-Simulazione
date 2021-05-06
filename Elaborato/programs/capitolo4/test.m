n = 100;
A = rand(n);
B = A+A';
C = eps*(A-A');
[D,~] = qr(A);

f = @(z) sin(z);
F1 = funm(C,@sin);
F2 = funm_normal(C,f);
e = norm(F1-F2)
