function [alpha,beta,order] = adams_bashforth_coefficients(k)
%ADAMS_BASHFORTH_COEFFICIENTS Calcola i coefficienti alpha e beta del metodo
% di Adams-Bashforth a k passi (il cui ordine è k).
    A = zeros(2*k+2,2*k+2);
    b = zeros(2*k+2,1);
    n = 1;
    % imponi alpha_k = 1
    A(n,k+1) = 1;
    b(n) = 1;
    n = n+1;
    % imponi alpha_{k-1} = -1
    A(n,k) = -1;
    b(n) = 1;
    n = n+1;
    % imponi che gli altri alpha siano nulli
    for i = 1:k-1
        A(n,i) = 1;
        n = n+1;
    end
    % imponi beta_k = 0
    A(n,2*k+2) = 1;
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


