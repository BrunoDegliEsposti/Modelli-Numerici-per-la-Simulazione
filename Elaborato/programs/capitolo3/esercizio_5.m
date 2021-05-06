%% Adams-Bashforth
for k = 1:6
    A = zeros(2*k+2,2*k+2);
    b = zeros(2*k+2,1);
    n = 1;

    % alpha_k = 1
    A(n,k+1) = 1;
    b(n) = 1;
    n = n+1;
    
    % alpha_{k-1} = -1
    A(n,k) = -1;
    b(n) = 1;
    n = n+1;
    
    % annulla gli altri alpha
    for i = 1:k-1
        A(n,i) = 1;
        n = n+1;
    end
    
    % annulla beta_k
    A(n,2*k+2) = 1;
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
    disp(fliplr(beta));
end

%% Adams-Moulton
for k = 1:6
    A = zeros(2*k+2,2*k+2);
    b = zeros(2*k+2,1);
    n = 1;

    % alpha_k = 1
    A(n,k+1) = 1;
    b(n) = 1;
    n = n+1;
    
    % alpha_{k-1} = -1
    A(n,k) = -1;
    b(n) = 1;
    n = n+1;
    
    % annulla gli altri alpha
    for i = 1:k-1
        A(n,i) = 1;
        n = n+1;
    end
    
    % annulla il coefficiente di h^j
    for j = 1:k+1
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
    disp(fliplr(beta));
end
