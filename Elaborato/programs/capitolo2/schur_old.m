n = 14;
syms x;
pse = poly2sym([1 -0.9],x);
for i = 2:n
    pse = pse * poly2sym([1 -0.9],x);
end
pn1 = sym2poly(pse);
pn2 = poly(0.9*ones(1,n));
psi1 = poly2sym(pn1,x);
psi2 = poly2sym(pn2,x);
r1 = max(abs(roots(pn1)));
r2 = max(abs(roots(pn2)));
digits(1000);
r3 = max(abs(double(vpasolve(psi1==0,x))));
r4 = max(abs(double(vpasolve(psi2==0,x))));
[r1,r2,r3,r4,is_schur(pn1),is_schur(pn2)]

n = 30;
syms x;
pse = poly2sym([1 -1/n],x);
for i = 2:n-1
    pse = pse * poly2sym([1 -i/n],x);
end
pn = sym2poly(pse);
psi = poly2sym(pn,x);
r1 = roots(pn);
digits(1000);
r2 = double(vpasolve(psi==0,x));
[max(abs(r1)),max(abs(r2)),is_schur(pn)]
