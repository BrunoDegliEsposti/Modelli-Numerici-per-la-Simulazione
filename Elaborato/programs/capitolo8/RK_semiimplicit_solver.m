function [Y,errcode] = RK_semiimplicit_solver(odefun,t,Y0,h,A_ii,g,toll)
%RK_SEMIIMPLICIT_SOLVER Metodo di Newton stazionario per l'i-esima equazione
% non lineare F(Y) = 0 che dev'essere risolta in un metodo Runge-Kutta
% semiimplicito, con F(Y) = Y - h*A_ii*odefun(t,Y) - g.
% Il tempo t dev'essere calcolato come t_n+h*c_i.
% Y0 è il punto di partenza per il metodo iterativo.
% Il criterio di arresto utilizzato è misto con tolleranza toll.
% Se tale criterio viene soddisfatto in meno di 250 iterazioni,
% allora errcode=0, altrimenti errcode=1.
    m = length(Y0);
    J0 = numerical_jacobian(odefun,t,Y0);
    [L,U] = lu(eye(m)-h*A_ii*J0);
    Y = Y0;
    for n = 1:250
        F = Y - h*A_ii*odefun(t,Y) - g;
        dY = U\(L\F);
        Y = Y - dY;
        if all(abs(dY) <= toll + toll*abs(Y))
            errcode = 0;
            return
        end
    end
    errcode = 1;
end
