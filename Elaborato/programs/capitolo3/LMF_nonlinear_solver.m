function [y,errcode] = LMF_nonlinear_solver(odefun,t0,y0,h,beta_k,g,toll)
%LMF_NONLINEAR_SOLVER Metodo di Newton stazionario per l'equazione
% non lineare G(y) = 0 che dev'essere risolta a ogni passo di
% un metodo LMF implicito, con G(y) = y - h*beta_k*odefun(t0,y) - g.
% y0 è il punto di partenza per il metodo iterativo.
% Il criterio di arresto utilizzato è misto con tolleranza toll.
% Se tale criterio viene soddisfatto in meno di 250 iterazioni,
% allora errcode=0, altrimenti errcode=1.
    m = length(y0);
    J0 = numerical_jacobian(odefun,t0,y0);
    [L,U] = lu(eye(m)-h*beta_k*J0);
    y = y0;
    for n = 1:250
        G = y - h*beta_k*odefun(t0,y) - g;
        dy = U\(L\G);
        y = y - dy;
        if all(abs(dy) <= toll + toll*abs(y))
            errcode = 0;
            return
        end
    end
    errcode = 1;
end
