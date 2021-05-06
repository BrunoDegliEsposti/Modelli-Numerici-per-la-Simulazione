function [F, c] = funm_naive(A, fun)
%FUN_NAIVE Calcola f(A) diagonalizzando A. Fornisce anche il numero di
% condizionamento della matrice di cambio di base. Se A è normale,
% l'algoritmo è stabile e produce lo stesso risultato di funm.
% Se A non è diagonalizzabile, viene calcolata la media tra f(A+R) 
% e f(A-R), con R matrice random di norma sqrt(eps) circa.
    [n,m] = size(A);
    if n ~= m
        error("La matrice in ingresso non è quadrata");
    end
    [U,T] = schur(A, 'complex');
    if norm(triu(T,1)) < 1e-15
        for i = 1:n
            T(i,i) = fun(T(i,i));
        end
        F = U*T*U';
        c = 1;
        return
    end
    [V,D] = eig(T);
    if rank(V) == n       
        for i = 1:n
            D(i,i) = fun(D(i,i));
        end
        F = U*V*D*inv(V)*U';
        c = cond(V);
    else
        R = sqrt(eps) * rand(n);
        [F1, c1] = funm_diag(T+R, fun);
        [F2, c2] = funm_diag(T-R, fun);
        F = U * (0.5*(F1+F2)) * U';
        c = max(c1,c2);
    end
    if c > 1e8
        warning('Poche cifre significative nel risultato');
    end
end


