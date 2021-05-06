function [yT,tab] = simpson_rule(odefun,t0,T,h,y0)
%SIMPSON_RULE Metodo di Cavalieri-Simpson con passo h per il problema
% di Cauchy y'(t) = odefun(t,y(t)), y(t0) = y0, t0 <= t <= T.
% L'argomento in uscita yT è l'approssimazione numerica di y(T).
% L'argomento in uscita tab è facoltativo, è di tipo "table" e contiene
% tutte le iterazioni effettuate dal metodo.
    save_table = (nargout > 1);
    alpha = [-1,0,1];
    beta = [1/3,4/3,1/3];
    order = 4;
    if save_table
        [yT,tab] = LMF(odefun,t0,T,h,y0,alpha,beta,order);
    else
        yT = LMF(odefun,t0,T,h,y0,alpha,beta,order);
    end
end
