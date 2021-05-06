function [y1,ML] = kepler( t, y )
%
%   Usage:
%   y' = kepler( t, y )  or  [H,ML] = kepler( y ) 
%
%   Kepler problem. First order formulation.  y = [q1 q2 p1 p2]'.
%   ML contains the angular momentum and the second component of 
%   the Lenz vector.
%
%   See pages 131-134  of
%   L.Brugnano, F.Iavernaro.
%   Line Integral Methods for Conservative Problems.
%   Chapman and Hall/CRC - 2015
%   Series: Monographs and Research Notes in Mathematics.
%
if nargin==0
   error('see kepler0')
elseif nargin==1   % Hamiltonian
   y    = t;
   sroq = sqrt( y(:,1).^2+y(:,2).^2 );
   rop  = y(:,3).^2+y(:,4).^2;
   y1   = .5*rop-1./sroq;
   if nargout==2  % angular momentum and Lenz vector
      ML = [y(:,1).*y(:,4)-y(:,2).*y(:,3) ...
            y(:,2).*y(:,3).^2-y(:,1).*y(:,3).*y(:,4)-y(:,2)./sroq];
   end
else % vector field
   roq1 = ( y(:,1).^2+y(:,2).^2 ).^(3/2);
   y1   = [y(:,3:4) -y(:,1)./roq1 -y(:,2)./roq1];
end
return      
   
