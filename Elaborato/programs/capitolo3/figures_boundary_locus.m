%% Adams-Bashforth
fig1 = figure(1);
fig1.Units = 'centimeters';
fig1.InnerPosition(3:4) = [14 14];
hold on;
plot([-10,10],[0,0],'k--');
plot([0,0],[-10,10],'k--');
for k=1:6
    [alpha,beta] = adams_bashforth_coefficients(k);
    p(k) = plot_boundary_locus(alpha,beta);
end
axis([-2.5,0.5,-1.5,1.5]);
xlabel('x');
ylabel('y');
legend(p,{'$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','$k=6$'},...
    'Interpreter','latex','Location','northwest');
exportgraphics(fig1,'../../figures/capitolo3/adams-bashforth-boundary-loci.pdf',...
    'ContentType','vector','BackgroundColor','none');

%% Adams-Moulton
fig2 = figure(2);
fig2.Units = 'centimeters';
fig2.InnerPosition(3:4) = [14 14];
hold on;
plot([-10,10],[0,0],'k--');
plot([0,0],[-10,10],'k--');
for k=1:6
    [alpha,beta] = adams_moulton_coefficients(k);
    p(k) = plot_boundary_locus(alpha,beta);
end
axis([-7.5,0.5,-4,4]);
xlabel('x');
ylabel('y');
legend(p,{'$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','$k=6$'},...
    'Interpreter','latex','Location','northwest');
exportgraphics(fig2,'../../figures/capitolo3/adams-moulton-boundary-loci.pdf',...
    'ContentType','vector','BackgroundColor','none');

%% BDF
fig3 = figure(3);
fig3.Units = 'centimeters';
fig3.InnerPosition(3:4) = [14 14];
hold on;
plot([-100,100],[0,0],'k--');
plot([0,0],[-100,100],'k--');
for k=1:6
    [alpha,beta] = BDF_coefficients(k);
    p(k) = plot_boundary_locus(alpha,beta);
end
axis([-10,30,-25,25]);
xlabel('x');
ylabel('y');
legend(p,{'$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','$k=6$'},...
    'Interpreter','latex','Location','northeast');
exportgraphics(fig3,'../../figures/capitolo3/BDF-boundary-loci.pdf',...
    'ContentType','vector','BackgroundColor','none');
