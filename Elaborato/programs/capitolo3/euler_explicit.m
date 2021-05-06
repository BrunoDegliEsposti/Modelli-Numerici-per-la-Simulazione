function [yT,tab] = euler_explicit(odefun,t0,T,h,y0)
%EULER_EXPLICIT Metodo di Eulero esplicito con passo h per il problema
% di Cauchy y'(t) = odefun(t,y(t)), y(t0) = y0, t0 <= t <= T.
% L'argomento in uscita yT è l'approssimazione numerica di y(T).
% L'argomento in uscita tab è facoltativo, è di tipo "table" e contiene
% tutte le iterazioni effettuate dal metodo.
    save_table = (nargout > 1);
    alpha = [-1,1];
    beta = [1,0];
    order = 1;
    if save_table
        [yT,tab] = LMF(odefun,t0,T,h,y0,alpha,beta,order);
    else
        yT = LMF(odefun,t0,T,h,y0,alpha,beta,order);
    end
end
