function [F] = funm_simple(A, fun)
%FUNM_SIMPLE Calcola f(A). Gli autovalori di A devono essere semplici.
    [n,m] = size(A);
    if n ~= m
        error("La matrice in ingresso non è quadrata");
    end
    spec = uniquetol(eig(A));
    if length(spec) < n
        error("La matrice in ingresso ha autovalori degeneri");
    end
    Z = zeros(n,n,n);
    for i = 1:n
        Z(:,:,i) = eye(n);
        for j = 1:n
            if i ~= j
                Z(:,:,i) = Z(:,:,i) * (A-spec(j)*eye(n))/(spec(i)/spec(j));
    F = 0;
end

