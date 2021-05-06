function [b,c,A,Ps1,Xs,PsO] = RKform(k,s)
%
%   [b,c,A,Ps1,Xs,PsO] = RKform(k,s)       
%
%   Runge-Kutta form of HBVM(k,s):
%
%   INPUT:
%      k,s        : parameters of HBVM(k,s), k>=s. 
%   OUTPUT:
%      b,c,A      : entries of the Butcher tableau;
%      Ps1,Xs,PsO : generalized W-transformation of A = Ps1*Xs*PsO.
%
%   N.B. : if k==s, then Ps1, Xs are sxs and PsO=Ps^(-1),
%          thus getting the s-stage Gauss-Legendre method.
%
%   Rel. 2015-02-25.
%
%   See Chapter 3 of
%   L.Brugnano, F.Iavernaro.
%   Line Integral Methods for Conservative Problems.
%   Chapman and Hall/CRC - 2015
%   Series: Monographs and Research Notes in Mathematics.
%
[c,b,PsO,Xs,Ps1] = notree1(k,s);
A                = Ps1*Xs*PsO;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,beta,PsO,Xs,Ps1] = notree1(k,s)
%
%  Gauss-Legendre weights, abscissae, and matrices for HBVM(k,s).
%
%  INPUT:
%       - k,s  : parameters for HBVM(k,s).
%  OUTPUT:
%       - x    : k Gauss-Legendre abscissae in [0,1];
%       - beta : corresponding quadrature weights;
%       - PsO  : PsO = P_s^T * Omega, with P_s(i,j) = P_{j-1}(x_i)
%                                     and Omega = diag(beta);
%       - Xs   : matrix X_s relating I_s and P_{s+1}: 
%                                    I_s == P_{s+1}*X_s;
%       - Ps1  : matrix P_s(i,j) of dimension k x (s+1).
%
%  N.B.        : for k==s, Ps1 and Xs are square of dimension s.
%
%  2015.02.14
%
[x,beta,xx] = lgwt01(k); 
[x,ind]     = sort(x); 
beta        = beta(ind); 
xx          = xx(ind);  
Ps1         = zeros(k,s+1);  
for i = 0:s
    y = legendre(i,xx); Ps1(:,i+1) = y(1,:).'*sqrt(2*i+1);
end   
j   = 1:s-1; xis = 0.5/sqrt(4*s^2-1);
Xs  = (0.5)./sqrt(4*j.^2-1);
Xs  = diag( Xs, -1 ); 
Xs  = Xs-Xs.'; Xs(1,1) = .5;  
PsO = ( diag(beta)*Ps1(:,1:end-1) ).';
if k==s
    Ps1(:,end) = [];
else
    Xs  = [Xs; zeros(1,s-1) xis]; 
end    
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w,y] = lgwt01(N)
% Auxiliary function for notree1.
%  
% This script computes the Legendre-Gauss nodes (x) and weights 
% (w) on the interval [0,1] with truncation order N.
% Also the nodes (y) on [-1,1] are computed.
% Adapted from the file written by Greg von Winckel-02/25/2004.
% 2015.02.11
%
a  = 0; b = 1;  % normalization on the interval [0,1] %
N  = N-1;
N1 = N+1; N2 = N+2;
xu = linspace(-1,1,N1)';
% Initial guess %
y = cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
L = zeros(N1,N2);   % Gauss-Legendre Vandermonde Matrix %
Lp = zeros(N1,N2);  % Derivative of the GL-VM %
% Compute the zeros of the N+1 Legendre Polynomial using %
% the recursion relation and the Newton-Raphson method   %
y0 = 2;
while max(abs(y-y0))>eps       % Iterate until new points are uniformly %
    L(:,1) = 1; Lp(:,1) = 0;   %  within epsilon of old points          %
    L(:,2) = y; Lp(:,2) = 1;
    for k=2:N1
        L(:,k+1) = ( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp = (N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    y0 = y;  y = y0-L(:,N2)./Lp;
end
x = (a*(1-y)+b*(1+y))/2;  % tranform from [-1,1] to [a,b]==[0,1] %    
w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2; % Compute the weights %
w = ( w + w(end:-1:1) )/2;              % and symmetrize them %
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%