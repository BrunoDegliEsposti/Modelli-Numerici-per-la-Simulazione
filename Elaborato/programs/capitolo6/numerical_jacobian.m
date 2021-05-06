function [Jf,f0] = numerical_jacobian(f,t0,y0)
%NUMERICAL_JACOBIAN Valuta la matrice jacobiana di f(t,y): R x R^m -> R^m
% rispetto all'argomento vettoriale y nel punto (t0,y0).
% La matrice Jf ha quindi dimensione (m,m). Il valore f0 è uguale a f(t0,y0).
% La matrice Jf viene calcolata tramite differenze finite del primo ordine,
% quindi ha errore relativo dell'ordine di sqrt(eps),
% a patto che le derivate del primo e secondo ordine in (t0,y0)
% siano trascurabili nell'analisi dell'errore.
    m = length(y0);
    Jf = zeros(m,m);
    f0 = f(t0,y0);
    h = sqrt(eps*sum(abs(y0)) + eps*sum(abs(f0)));
    h = max(1e-12,h);
    for i = 1:m
        y1 = y0;
        y1(i) = y1(i) + h;
        f1 = f(t0,y1);
        Jf(:,i) = (f1-f0)/h;
    end
end
