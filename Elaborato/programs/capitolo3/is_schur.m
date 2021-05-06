function [b] = is_schur(p)
% IS_SCHUR Primo criterio di Schur per la stabilità di un polinomio.
% p: polinomio a coefficienti complessi.
% b: booleano che indica se p ha solo radici di modulo < 1.
    [m,n] = size(p);
    assert(m==1);
    if n == 1
        b = true;
        return
    end
    q = conj(fliplr(p));
    p0 = p(end);
    q0 = q(end);
    if abs(q0) <= abs(p0)
        b = false;
        return
    end
    reduced = q0*p(1:end-1)-p0*q(1:end-1);
    reduced = reduced/median(abs(reduced));
    b = is_schur(reduced);
end
