S = zeros(55);
S(:,1) = 1;
n = 2;
exact = true;
while exact
    for i = 2:n
        S(n,i) = S(n-1,i-1) + i*S(n-1,i);
        if S(n,i) > 2^53
            exact = false;
        end
    end
    n = n+1;
end
S = S(1:n-2, 1:n-2);

% Il più grande numero di Stirling
% di seconda specie calcolato è S(22,9)