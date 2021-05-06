%% Paziente non diabetico
t1 = [0,   30, 60, 90, 120, 180];
G1 = [72, 130, 95, 52,  49,  74];

type1 = fittype('Gbar + A*exp(-alpha*t)*sin(omega*t)',...
    'dependent',{'G'},'independent',{'t'},...
    'coefficients',{'Gbar','A','alpha','omega'});
fit1 = fit(t1',G1',type1,'StartPoint',[70,60,0.01,2*pi/150]);

fig1 = figure(1);
fig1.Units = 'centimeters';
fig1.InnerPosition(3:4) = [18 12];
plot(t1,G1,'k+','MarkerSize',12);
hold on;
plot(fit1,'b');
xlabel('t (minuti)');
ylabel('G(t) (mg/dl)');
legend('Misure di G(t)','Fit ai minimi quadrati');
exportgraphics(fig1,'../../figures/capitolo6/diabete1.pdf',...
    'ContentType','vector','BackgroundColor','none');


%% Paziente diabetico
t2 = [0,    30,  60,  90, 120, 150, 180, 210, 240];
G2 = [136, 160, 195, 190, 180, 170, 160, 150, 140];

type2 = fittype('Gbar + A*exp(-alpha*t)*sin(omega*t)',...
    'dependent',{'G'},'independent',{'t'},...
    'coefficients',{'Gbar','A','alpha','omega'});
fit2 = fit(t2',G2',type2,'StartPoint',[140,60,0.01,2*pi/300]);

fig2 = figure(2);
fig2.Units = 'centimeters';
fig2.InnerPosition(3:4) = [18 12];
plot(t2,G2,'k+','MarkerSize',12);
hold on;
plot(fit2,'b');
xlabel('t (minuti)');
ylabel('G(t) (mg/dl)');
legend('Misure di G(t)','Fit ai minimi quadrati');
exportgraphics(fig2,'../../figures/capitolo6/diabete2.pdf',...
    'ContentType','vector','BackgroundColor','none');
