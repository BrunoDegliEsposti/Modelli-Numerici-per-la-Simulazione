%% Grafico

z = roots(poly(0.9*ones(1,13)));
plot(nsidedpoly(1024));
hold on;
scatter(real(z),imag(z),'xk');
axis equal;
axis([0.65,1.15,-0.25,0.25]);
ax = gca;
exportgraphics(ax,'../../figures/capitolo2/qr13.pdf',...
    'ContentType','vector','BackgroundColor','none');

%% Programma

% Definizione esatta del polinomio
k = 14;
syms z;
p1 = (z-0.9)^k;
% Arrotondamento dei coefficienti
p2 = sym2poly(p1);
% Calcolo della radice di modulo massimo
digits(1000);
p3 = poly2sym(p2,z);
zmax = max(abs(double(vpasolve(p3==0,z))));
