% definizione dei parametri e caricamento dei dati
tau = 5;
L = 120;
m = L/tau;
T = readtable('popolazione_italiana_2019.xls');
writetable(T,'../../tables/capitolo5/popolazione_italiana_2019.dat');
alpha = T.alpha;
beta = T.beta(1:m-1);
A = diag(beta,-1);
A(1,:) = alpha;

% sostituzioni in avanti
nmax = 16;
Y = zeros(m,nmax+1);
Y(:,1) = T.y;
for n = 1:nmax
    Y(:,n+1) = A*Y(:,n);
end

%% visualizzazione dei risultati
fig1 = figure(1);
fig1.Units = 'centimeters';
fig1.InnerPosition(3:4) = [30 30];
x = 0:5:115;
bar3(x,Y);
years = 19:5:99;
xticklabels(num2str(years'));
xlabel('Anno (20xx)'); xticks(1:17);
ylabel('Età'); yticks(x); ylim([-5,120]);
zlabel('Individui');
% exportgraphics(fig1,'../../figures/capitolo5/leslie-dinamica-3D.png',...
%     'ContentType','image','Resolution',300,'BackgroundColor','white');

%% autovalore e autovettore dominanti
fig2 = figure(2);
fig2.Units = 'centimeters';
fig2.InnerPosition(3:4) = [15 12];
[v0,lambda0] = eigs(A,1);
v0 = v0 ./ sum(v0);
Ynormalized = Y ./ sum(Y);

hold on;
p1 = plot(x,100*Ynormalized(:,1),'Color',[0,0,1],'LineWidth',1);
for n = 3:2:nmax-1
    rgb = [(n-1)/nmax,0,1-(n-1)/nmax];
    plot(x,100*Ynormalized(:,n),'Color',rgb,'LineWidth',1);
end
p2 = plot(x,100*Ynormalized(:,nmax+1),'Color',[1,0,0],'LineWidth',1);
p3 = plot(x,100*v0,'k','LineWidth',2);

xlabel("Fascia di et\`{a}",'Interpreter','latex'); xticks(x); xtickangle(45); xlim([0,105]);
ylabel('Percentuale di individui','Interpreter','latex');
grid on;
leg = legend([p1,p2,p3],...
    {'Popolazione iniziale $y(0)$','Popolazione finale $y(16)$','Autovettore dominante $v_0$'},...
    'Interpreter','latex');
leg.Units = 'centimeters';
leg.Position(1:2) = [2.85 2];
exportgraphics(fig2,'../../figures/capitolo5/leslie-dinamica-normalizzata.pdf',...
    'ContentType','vector','BackgroundColor','none');




