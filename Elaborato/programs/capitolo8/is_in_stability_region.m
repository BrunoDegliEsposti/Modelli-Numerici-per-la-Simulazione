function [is_inside] = is_in_stability_region(q,A,b)
%IS_IN_STABILITY_REGION Controlla se i punti q appartengono alla
% regione di assoluta stabilità di un metodo Runge-Kutta con pesi b e
% matrice di Butcher A.
    [m,n] = size(q);
    s = size(A,1);
    e = ones(s,1);
    is_inside = zeros(m,n,'logical');
    for j = 1:n
        for i = 1:m
            d1 = det(eye(s)-q(i,j)*A+q(i,j)*e*b(:)');
            d2 = det(eye(s)-q(i,j)*A);
            is_inside(i,j) = (abs(d1) < abs(d2));
        end
    end
end
