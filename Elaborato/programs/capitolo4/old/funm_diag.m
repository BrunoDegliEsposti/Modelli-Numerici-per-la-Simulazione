function [F, c] = funm_diag(A, fun)
%FUNM_DIAG Funzione ausiliaria che calcola f(A) diagonalizzando A.
    n = length(A);
    [V,D] = eig(A);
    for i = 1:n
        D(i,i) = fun(D(i,i));
    end
    F = V*D*inv(V);
    c = cond(V);
end

