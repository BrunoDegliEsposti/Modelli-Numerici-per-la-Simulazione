A = [0,  0,  0,  0;
     0.5,0,  0,  0;
     0,  0.5,0,  0;
     0,  0,  1,  0];
b = [1/6,1/3,1/3,1/6];
x = linspace(-4,2,500);
y = linspace(-4,4,500);
[X,Y] = meshgrid(x,y);
Z = is_in_stability_region(X+1i*Y,A,b);
contourf(X,Y,Z,[1,1]);

hold on
plot([-4,2],[0,0],'k--');
plot([0,0],[-4,4],'k--');
xlabel('x'); xticks(-4:0.5:2);
ylabel('y');
fig1 = gcf;
fig1.Units = 'centimeters';
fig1.InnerPosition(3:4) = [15 10];
exportgraphics(fig1,'../../figures/capitolo8/regione-RK-classico.png',...
    'ContentType','image','Resolution',300,'BackgroundColor','white');