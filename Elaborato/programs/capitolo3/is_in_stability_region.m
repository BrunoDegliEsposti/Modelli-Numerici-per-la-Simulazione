function [is_inside] = is_in_stability_region(q,alpha,beta)
%IS_IN_STABILITY_REGION Controlla se i punti q appartengono alla regione
% di assoluta stabilità di un metodo LMF con coefficienti alpha e beta.
    [m,n] = size(q);
    rho = fliplr(alpha);
    sigma = fliplr(beta);
    is_inside = zeros(m,n,'logical');
    for j = 1:n
        for i = 1:m
            is_inside(i,j) = is_schur(rho-q(i,j)*sigma);
        end
    end
end
