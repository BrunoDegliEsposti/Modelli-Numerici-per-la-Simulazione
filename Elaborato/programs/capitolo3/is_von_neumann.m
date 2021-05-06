function [b] = is_von_neumann(p)
% IS_VON_NEUMANN Secondo criterio di Schur per la stabilità di un polinomio.
% p: polinomio a coefficienti complessi in aritmetica ESATTA.
% b: booleano che indica se p è un polinomio di Von Neumann.
    [m,n] = size(p);
    assert(m==1 && isa(p,'sym'));
    if n == 1
        b = true;
        return
    end
    q = conj(fliplr(p));
    p0 = p(end);
    q0 = q(end);
    reduced = q0*p(1:end-1)-p0*q(1:end-1);
    dp = p(1:end-1).*(n-1:-1:1);
    b = (all(reduced==0) && is_schur(dp)) || ...
        (abs(q0) > abs(p0) && is_von_neumann(reduced));
    return
end
