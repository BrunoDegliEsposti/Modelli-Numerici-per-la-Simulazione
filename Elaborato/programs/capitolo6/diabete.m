% Paziente non diabetico
t1 = [0,   30, 60, 90, 120, 180];
G1 = [72, 130, 95, 52,  49,  74];
type1 = fittype('Gbar + A*exp(-alpha*t)*sin(omega*t)',...
    'dependent',{'G'},'independent',{'t'},...
    'coefficients',{'Gbar','A','alpha','omega'});
fit1 = fit(t1',G1',type1,'StartPoint',[70,60,0.01,2*pi/150]);
Gbar1 = fit1.Gbar; T1 = 2*pi/fit1.omega/60;

% Paziente diabetico
t2 = [0,    30,  60,  90, 120, 150, 180, 210, 240];
G2 = [136, 160, 195, 190, 180, 170, 160, 150, 140];
type2 = fittype('Gbar + A*exp(-alpha*t)*sin(omega*t)',...
    'dependent',{'G'},'independent',{'t'},...
    'coefficients',{'Gbar','A','alpha','omega'});
fit2 = fit(t2',G2',type2,'StartPoint',[140,60,0.01,2*pi/300]);
Gbar2 = fit2.Gbar; T2 = 2*pi/fit2.omega/60;
