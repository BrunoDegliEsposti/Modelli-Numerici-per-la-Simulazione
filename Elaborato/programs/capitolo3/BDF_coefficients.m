function [alpha,beta,order] = BDF_coefficients(k)
%BDF_COEFFICIENTS Calcola i coefficienti alpha e beta della backward
% differentiation formula a k passi (il cui ordine è k).
    A = zeros(2*k+2,2*k+2);
    b = zeros(2*k+2,1);
    n = 1;
    % imponi beta_k = 1
    A(n,2*k+2) = 1;
    b(n) = 1;
    n = n+1;
    % imponi che gli altri beta siano nulli
    for i = k+2:2*k+1
        A(n,i) = 1;
        n = n+1;
    end
    % imponi che la somma degli alpha sia nulla
    A(n,1:k+1) = 1;
    n = n+1;
    % imponi che il metodo abbia ordine k
    for j = 1:k
        for i = 0:k
            A(n,i+1) = i^j;
            A(n,k+1+i+1) = -j*(i^(j-1));
        end
        n = n+1;
    end
    % calcola x in aritmetica esatta
    x = sym(A)\sym(b);
    alpha = double(x(1:k+1)');
    beta = double(x(k+2:2*k+2)');
    order = k;
end
