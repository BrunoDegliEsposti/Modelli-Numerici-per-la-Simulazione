function [v] = dominant_eigenvalue(A)
%DOMINANT_EIGENVALUE Approssimazione dell'autovalore dominante di A.
    [m,~] = size(A);
    v = rand(m,1);
    for i = 1:50
        v = A*v;
        v = v/norm(v);
    end
end
