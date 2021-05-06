a = 4;
y0 = 1/3;
nmax = 1e6;

y = zeros(nmax,1);
y(1) = y0;
for i = 2:nmax
    y(i) = a*y(i-1)*(1-y(i-1));
end

nbins = 100;
idx = discretize(y,nbins);
density = zeros(nbins,1);
for i = 1:nmax
    density(idx(i)) = density(idx(i)) + 1;
end
density = nbins*density/nmax;

t = linspace(0,1,nbins+1);
density_exact = 1./(pi*sqrt(t.*(1-t)));

fig = figure(1);
fig.Units = 'centimeters';
fig.InnerPosition(3:4) = [18 12];
plot(t,density_exact);
hold on;
stairs(t,[density;density(end)]);
xlabel('y');
ylabel('Densità');
legend('Densità esatta h(y)','Densità numerica (1e6 iterazioni)',...
    'Location','southeast');
ylim([0,3]);
exportgraphics(fig,'../../figures/capitolo7/logistic-ergodicity.pdf',...
    'ContentType','vector','BackgroundColor','none');



