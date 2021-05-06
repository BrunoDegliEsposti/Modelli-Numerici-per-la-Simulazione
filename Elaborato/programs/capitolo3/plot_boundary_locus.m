function [p] = plot_boundary_locus(alpha,beta)
%PLOT_BOUNDARY_LOCUS Disegna il boundary locus del metodo LMF
% definito dai coefficienti alpha e beta.
    rho = fliplr(alpha);
    sigma = fliplr(beta);
    theta = linspace(0,2*pi,1000);
    z = exp(1i*theta);
    rhoz = polyval(rho,z);
    sigmaz = polyval(sigma,z);
    q = rhoz./sigmaz;
    mask = (sigmaz==0);
    q(mask) = NaN;
    p = plot(real(q),imag(q));
end
