function [t, y] = lorenz_solver(t_trans, t_len, y0, sigma, beta, rho)
    %LORENZ_SOLVER Simula il sistema dinamico di Lorenz
    %   Descrizione...
    f = @(t,y)[sigma*(y(2)-y(1));
               rho*y(1)-y(2)-y(1)*y(3);
               y(1)*y(2)-beta*y(3)];
    options = odeset('RelTol',1e-10,'AbsTol',1e-12);
    if t_trans ~= 0         % a ode45 non piace tspan == [0,0]
        [~, y] = ode45(f, [0 t_trans], y0, options);
        [t, y] = ode45(f, [0 t_len], y(end,:), options);
    else
        [t, y] = ode45(f, [0 t_len], y0, options);
    end
end

