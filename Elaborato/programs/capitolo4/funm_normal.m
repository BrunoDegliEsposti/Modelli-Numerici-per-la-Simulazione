function [fA] = funm_normal(A,f)
%FUNM_NORMAL Calcola la funzione di matrice f(A) nel caso in cui A sia una
% matrice normale (A commuta con la sua trasposta coniugata).
% L'argomento in ingresso f:C->C è un function handle.
    [mbar,ncols] = size(A);
    if mbar ~= ncols
        error('La matrice in ingresso non è quadrata');
    end
    [Q,T] = schur(A,'complex');
    D = diag(T);
    nd = max(abs(D));
    nu = max(max(abs(triu(T,1))));
    if nu > 1e-14 * nd
        error('La matrice in ingresso non è normale');
    end
    for i=1:mbar
        D(i) = f(D(i));
    end
    fA = Q*diag(D)*Q';
end
