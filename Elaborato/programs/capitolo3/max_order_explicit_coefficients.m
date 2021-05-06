function [alpha_exact,beta_exact] = max_order_explicit_coefficients(k)
%MAX_ORDER_EXPLICIT_COEFFICIENTS Calcola in aritmetica esatta i coefficienti
% alpha e beta dell'unico metodo lineare a k passi esplicito di ordine 2k-1.
% A e b definiscono un sistema lineare di vincoli la cui soluzione x = A\b
% contiene i coefficienti alpha_0, ..., alpha_k, beta_0, ..., beta_k
% con indici                 1,    ...,    k+1,    k+2,  ..., 2*k+2.
% La variabile n tiene traccia del numero del vincolo che si va a imporre.
    A = zeros(2*k+2,2*k+2);
    b = zeros(2*k+2,1);
    n = 1;
    % imponi alpha_k = 1
    A(n,k+1) = 1;
    b(n) = 1;
    n = n+1;
    % imponi beta_k = 0
    A(n,2*k+2) = 1;
    n = n+1;
    % imponi che la somma degli alpha sia nulla
    A(n,1:k+1) = 1;
    n = n+1;
    % imponi che il metodo abbia ordine 2k-1
    for j = 1:2*k-1
        for i = 0:k
            A(n,i+1) = i^j;
            A(n,k+1+i+1) = -j*(i^(j-1));
        end
        n = n+1;
    end
    % calcola alpha e beta in aritmetica esatta
    x = sym(A)\sym(b);
    alpha_exact = x(1:k+1)';
    beta_exact = x(k+2:2*k+2)';
end
