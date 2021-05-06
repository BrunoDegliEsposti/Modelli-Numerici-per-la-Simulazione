function [yT,tab] = adams_bashforth(odefun,t0,T,h,y0,k)
%ADAMS_BASHFORTH Metodo di Adams-Bashforth a k passi e passo temporale h per
% il problema di Cauchy y'(t) = odefun(t,y(t)), y(t0) = y0, t0 <= t <= T.
% L'argomento in uscita yT è l'approssimazione numerica di y(T).
% L'argomento in uscita tab è facoltativo, è di tipo "table" e contiene
% tutte le iterazioni effettuate dal metodo.
    save_table = (nargout > 1);
    [alpha,beta,order] = adams_bashforth_coefficients(k);
    if save_table
        [yT,tab] = LMF(odefun,t0,T,h,y0,alpha,beta,order);
    else
        yT = LMF(odefun,t0,T,h,y0,alpha,beta,order);
    end
end
