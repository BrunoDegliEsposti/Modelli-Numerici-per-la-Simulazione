function [sigma,kappa,gamma] = stiffness_ratio(f,t0,T,y0)
%STIFFNESS_RATIO Stima del rapporto di stiffness per il problema
% differenziale y'(t) = f(t,y(t)) con t in [t0,T] e y(t0) = y0.
    
    % stima della direzione ottimale lungo cui perturbare y0
    J = numerical_jacobian(f,t0,y0);
    v = dominant_eigenvalue(J);
    y0_perturbed = y0 + 1e-8*v;
    
    % soluzione "esatta" e perturbata
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    sol1 = ode15s(f,[t0,T],y0,opts);
    sol2 = ode15s(f,[t0,T],y0_perturbed,opts);
    
    % interpolazione delle soluzioni su un dominio comune
    t = sort(unique([sol1.x,sol2.x]));
    y1 = deval(sol1,t);
    y2 = deval(sol2,t);    
    
    % stima del rapporto di stiffness
    dy = vecnorm(y1-y2);
    kappa = max(dy);
    gamma = trapz(t,dy)/(T-t0);
    sigma = kappa/gamma;
end
