function [Jf,f0] = numjacobian(fun,t0,y0)
%
%  NUMERICAL JACOBIAN  (slightly adapted from the code BiM/BiMD).
%  2015.02.12
%
    m  = length(y0);
    Jf = zeros(m);
    f0 = reshape( feval( fun, t0, y0 ), m, 1 );
    dd = eps^(1/3);
    for i = 1:m
        ysafe = y0(i);
        delt  = sqrt( eps*max(dd,abs(ysafe)) );
        y0(i) = ysafe +delt;
        f1    = feval( fun, t0, y0 );
        Jf(:,i) = ( f1(:) -f0 )/delt;
        y0(i) = ysafe;
    end
    if nargout>1
        f0 = f0.';
    end
end
