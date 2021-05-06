%% BDF
for k = 1:6
    A = zeros(2*k+2,2*k+2);
    b = zeros(2*k+2,1);
    n = 1;
    
    % beta_k = 1
    A(n,2*k+2) = 1;
    b(n) = 1;
    n = n+1;
    
    % annulla gli altri beta
    for i = k+2:2*k+1
        A(n,i) = 1;
        n = n+1;
    end
    
    % annulla la somma degli alpha
    A(n,1:k+1) = 1;
    n = n+1;
    
    % annulla il coefficiente di h^j
    for j = 1:k
        for i = 0:k
            A(n,i+1) = i^j;
            A(n,k+1+i+1) = -j*(i^(j-1));
        end
        n = n+1;
    end

    % calcola x in aritmetica esatta
    x = sym(A)\sym(b);
    alpha = x(1:k+1)';
    beta = x(k+2:2*k+2)';
    rho = fliplr(alpha);
    fprintf('La BDF di grado %d\n',k);
    fprintf('rho = %s\n', char(rho));
    if is_von_neumann(rho)
        fprintf('è 0-stabile\n\n');
    else
        fprintf('non è 0-stabile\n\n');
    end
end







