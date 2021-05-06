function [Y,errcode] = RK_implicit_solver(odefun,t,y,Y0,h,A,c,toll)
%RK_IMPLICIT_SOLVER Metodo di Newton stazionario per l'equazione
% non lineare F(Y) = 0 che dev'essere risolta in un metodo Runge-Kutta
% implicito, con F(Y) = Y - e \otimes y - h*(A \otimes I_m)*fY.
% Il vettore fY è definito come l'incolonnamento di odefun(t+h*c_i,Y_i)
% per i che va da 1 a s. Y0 è il punto di partenza per il metodo iterativo.
% Il criterio di arresto utilizzato è misto con tolleranza toll.
% Se tale criterio viene soddisfatto in meno di 250 iterazioni,
% allora errcode=0, altrimenti errcode=1.
    s = size(A,1);
    m = length(y);
    J0 = numerical_jacobian(odefun,t,y);
    [L,U] = lu(eye(s*m)-h*kron(A,J0));
    Y = Y0;
    fY = zeros(m,s);
    A_otimes_Im = kron(A,eye(m));
    e_otimes_y = kron(ones(s,1),y);
    for n = 1:250
        for i = 1:s
            fY(:,i) = odefun(t+h*c(i),Y((1:m)+m*(i-1)));
        end
        F = Y - h*A_otimes_Im*fY(:) - e_otimes_y;
        dY = U\(L\F);
        Y = Y - dY;
        if all(abs(dY) <= toll + toll*abs(Y))
            errcode = 0;
            return
        end
    end
    errcode = 1;
end
