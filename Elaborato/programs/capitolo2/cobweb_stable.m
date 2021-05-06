%% Calcola una dinamica stabile

d0 = 10;
s0 = 0;
a = 1.5;
b = 1;
nmax = 10;
p = zeros(nmax,1);
p(1) = 3.2;
for n = 2:nmax
    p(n) = (b*p(n-1)+s0-d0)/(-a);
end

%% Disegna il grafico

fig = figure(1);
fig.Units = 'centimeters';
fig.InnerPosition(3:4) = [10 8];
hold on;
x = linspace(3,5,100);
supply = plot(x,b*x+s0,'b','LineWidth',1);
demand = plot(x,-a*x+d0,'r','LineWidth',1);
for n = 2:nmax
    x = p(n-1:n);
    y = [b*p(n-1)+s0, -a*p(n)+d0];
    %plot(x,y,'--k');
    quiver(x(1),y(1),0.95*(x(2)-x(1)),0.95*(y(2)-y(1)),0,'k');
end
for n = 2:nmax
    x = [p(n),p(n)];
    y = [-a*p(n)+d0,b*p(n)+s0];
    %plot(x,y,'--k');
    quiver(x(1),y(1),0.95*(x(2)-x(1)),0.95*(y(2)-y(1)),0,'k');
end
hold off;
%axis padded;
xlabel('Prezzo $p$','Interpreter','latex');
leg = legend([supply,demand],{'Offerta $S(p)$','Domanda $D(p)$'},'Interpreter','latex');
leg.Units = 'centimeters';
leg.Position(1:2) = [5.25 6];
exportgraphics(fig,'../../figures/capitolo2/cobweb1.pdf',...
    'ContentType','vector','BackgroundColor','none');





