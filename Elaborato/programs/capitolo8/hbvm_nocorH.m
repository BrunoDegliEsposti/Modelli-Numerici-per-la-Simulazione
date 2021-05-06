function [y,t,itot,ierr,etim]=hbvm(sflag, fun, y0, k, s, h, n, nob, stride)
%__________________________________________________________________________
%
% [y,t,it,ierr,etim] = hbvm( sflag, fun, y0, k, s, h, n, nob, stride )
%
%               Fixed stepsize implementation of the  HBVM(k,s)  method for 
%               canonical Hamiltonian problems with autonomous Hamiltonian.
%               Input problem either in first order form:
%
%                   q' = nabla_p H(q,p),     p' = -nabla_q H(q,p),      (1)
%
%               with Hamiltonian  H(q,p),  or  second order form:
%
%                                 q'' = M * nabla U(q),                 (2)
%
%               with Hamiltonian  H(q,p) = 1/2 p^T*M*p -U(q)  and  M = M^T.
%
%               The code may also handle general ODE problems, 
%                   y' = f(t,y)    or    y'' = f(t,y)
%               provided that fun is set as below specified.
%
%  INPUT PARAMETERS:
%
%   sflag = flag for first order problems (0) or second order problems 
%                (if different from 0).
%
%   fun   = string containing the idetifier of the function of the problem
%           which must be a canonical (autonomous) Hamiltonian one, i.e.,
%
%           For the first order problem (1):
%            fun( t, [q,p] ) -> returns the rhs of the ODE
%                               (f(t,[q,p]), for non Hamiltonian problems);
%            fun( [q,p] )    -> returns the Hamiltonian 
%                               (void, [], for non Hamiltonian problems).
%           For the second order problem (2):
%            fun( t, q )     -> returns nabla U(q); 
%                               (f(t,q), for non Hamiltonian problems);
%            fun()           -> returns M
%                               (M = 1, for non Hamiltonian problems);
%            fun( [q,p] )    -> returns the Hamiltonian 
%                               (void, [], for non Hamiltonian problems).
%
%    N.B.: fun(t,[q,p]) / fun(t,q) has to work also with vector inputs, and
%    ====> has to return a row for each evaluation of the rhs of the ODE.
%
%   y0    = initial condition == [q0,p0]. 
%
%   (k,s) = parameters defining the HBVM(k,s) method, k>=s. 
%           When k==s, one obtains the s-stage Gauss-Legendre RK method.
%
%   h     = time step.
%
%   n     = number of integration steps.
%
%   nob   = optional parameter which may be a number or a string:
%         if a number, the discrete problem is solved via a fixed-point
%                      iteration;
%         if a string, the discrete problem is solved via a blended
%                      iteration, with analytical jacobian, whose function
%                      is the specified string. Such a function has to be
%                      in the form: 
%                  G = jaco( t, [q,p] );   in case of first order problems;
%                  G = nabla^2 U( t, q ); in case of second order problems;
%         if void or not specified, the discrete problem is solved via a 
%                      blended iteration with numerical Jacobian (default).
%
%   stride = stride for the output of the computed solution. If stride is 
%            specified and positive, n/stride should be an integer.
%            If stride<=0, only the final point of the computed solution 
%               is provided in output.
%            The default value is stride=1.
%
% OUTPUT PARAMETERS:
%
%   (t,y) = computed solution (t is optional), spaced every stride steps,
%           y = [q,p].
%
%   itot  = total number of nonlinear iterations for solving the discrete 
%           problems.
%
%   ierr  = error flag: 
%           ==  0    if OK; 
%           ==  nu>0 if problems occurred with the nonlinear iteration at 
%                    nu integrations steps;
%           == -1    if NaNs are generated. In such a case, the integration
%                    is stopped at the current point.
%
%   etim  = execution time.
%
%                                                          Rel: 2020-02-13.
%__________________________________________________________________________
%                                                          
% References:
% -----------
% 1) L.Brugnano, F.Iavernaro, D.Trigiante. 
%    Hamiltonian Boundary Value Methods (Energy Preserving Discrete 
%    Line Integral Methods). 
%    Journal of Numerical Analysis, Industrial and Applied Mathematics 
%    5,1-2 (2010) 17-37.
% 2) L.Brugnano, F.Iavernaro, D.Trigiante. 
%    A note on the efficient implementation of Hamiltonian BVMs. 
%    Journal of Computational and Applied Mathematics 
%    236 (2011) 375-383.
% 3) L.Brugnano, F.Iavernaro, D.Trigiante.  
%    A simple framework for the derivation and analysis of effective 
%    one-step methods for ODEs. 
%    Applied Mathematics and Computation 218 (2012) 8475-8485.
% 4) L.Brugnano, F.Iavernaro.
%    Line Integral Methods for Conservative Problems.
%    Chapman and Hall/CRC - 2016
%    Series: Monographs and Research Notes in Mathematics.
% 5) L.Brugnano, F.Iavernaro.
%    Line Integral Solution of Differential Problems. 
%    Axioms 7(2) (2018) 36.
% 5) P.Amodio, L.Brugnano, F.Iavernaro.
%    Analysis of Spectral Hamiltonian Boundary Value Methods (SHBVMs) 
%    for the numerical solution of ODE problems.
%    Numerical Algorithms (2019) https://doi.org/10.1007/s11075-019-00733-7
%__________________________________________________________________________
%                                                          Rel: 2020-02-13.
% 
linea = '--------------------------------------------'; 
warning off backtrace
lib   = HBVMlib1;  % HBVMlib1 required % 
disp(linea);
fprintf('HBVM(%i,%i) \n',k,s);
disp(linea); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK INPUT
if nargin<7, error('not enough input arguments provided'); end
if ( k<s || s<=0 || k<=0 || n<=0 ), error('wrong input parameters'); end
corH = k>s;
if nargin<=8, stride = 1; end  % output at each step %
if nargin<=7, nob    = 0; end  % default, blended iteration, no jacobian %
if isempty(nob), nob = 0; end
if nob==0
   jaco = [];  nob =  0;  disp('blended iteration, numerical Jacobian');
elseif ischar(nob) 
   jaco = nob; nob = -1;  disp('blended iteration, analytical Jacobian');
else               
   jaco = [];  nob =  1;  disp('fixed-point iteration');
end    
if stride>0
   if fix(stride)<stride 
      error('stride has to be an integer')
   else 
      ns = fix(n/stride);
      if n/stride>ns 
         warning('the number of steps is not a multiple of stride') 
      end
   end
else
   stride = n; ns = 1;
end
if stride==1, disp('output at each step'); 
else,         fprintf('output every %i steps \n',stride);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if sflag==0   % first order problem %
   disp('first order problem');           disp(linea); 
   [y,t,itot,ierr,conta,steps,gs] = ... 
             lib.bhbvm1( fun, y0, k, s, h, n, nob, stride, jaco, ns, lib );
else          % second order problem %
   disp('special second order problem');  disp(linea); 
   [y,t,itot,ierr,conta,steps,gs] = ...
             lib.bhbvm2( fun, y0, k, s, h, n, nob, stride, jaco, ns, lib );
end
etim = toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK OUTPUT
if nob==1
   fprintf('number of fixed-point iterations  = %i \n',itot); 
else
   fprintf('number of blended iterations  = %i \n',itot);
end   
fprintf('elapsed time = %g sec \n',etim);
fprintf('|gamma_{s-1}| = %g \n',gs);
if corH && (conta>0)
   disp(linea) 
   fprintf('Hamiltonian error too large in \n      %i steps of %i. \n', ...
            conta, steps );
   fprintf('Consider to increase k. \n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(linea)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         HBVMlib1 (library for HBVMs)                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  h = HBVMlib1
%
%  Makes available the following functions for HBVMs:
%
%  - notree     % matrices and vectors needed for HBVM(k,s) methods %
%  - computeL   % approximate the Lipschitz constant and, % 
%               %   if required, the Jacobian transposed. %
%  - correggiH  % improves the Hamiltonian value, if within RIEs, %
%               %      for first order problems.                  %
%  - correggiHs % improves the Hamiltonian value, if within RIEs, %
%               %      for problems in the form q''=M*nabla U(q). %
%  - itera      % nonlinear iteration for first order problems. %
%  - iteras     % nonlinear iteration for problems %
%               %    in the form q''=M*nabla U(q). %
%  - bhbvm1     % HBVM(k,s) for first order problems %
%  - bhbvm2     % HBVM(k,s) for second order problems %
%
%                      / fun(t,y)  evaluates the rhs of the ODE;
%  Functions specifics:  fun(y)    evaluates the Hamiltonian;
%                      \ jaco(t,y) evaluates the Jacobian. If void,
%                                  this is numerically approximated.
%  For special second order problems in the form q'' = M*nabla U(q),
%                        fun(t,q)  evalates nabla U(q), and
%                        fun()     evaluates the (symmetric) matrix M;
%                        N.B.: if M = alfa * I, then fun may preferably
%                                               return the scalar alfa.
%
%  More details available in the headings of each function.
%
%  2015-04-18.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    h = HBVMlib1
h.notree      = @notree;     % matrices and vectors needed by HBVM(k,s) %
h.computeL    = @computeL;   % compute the Lip. constant and the Jacobian %
h.correggiH   = @correggiH;  % corrects the Hamiltonian, if within RIEs %
h.correggiHs  = @correggiHs; % same as above, for q''=M*nabla U(q) %
h.itera       = @itera;      % nonlinear iteration procedure %
h.iteras      = @iteras;     % same as above, for q''=M*nabla U(q) %
h.bhbvm1      = @bhbvm1;     % HBVM(k,s) for first order problems %
h.bhbvm2      = @bhbvm2;     % same as above, for q''=M*nabla U(q) %
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,t,itot,ierr,conta,steps,gs] = ...
                  bhbvm1( fun, y0, k, s, h, n, nob, stride, jaco, ns, lib )
%__________________________________________________________________________
%
%  Fixed stepsize implementation of the  HBVM(k,s) method for first order
%  canonical Hamiltonian problems with autonomous Hamiltonian.
%  nob =  0: blended iteration, numerical Jacobian;
%  nob = -1: blended iteration, analytical Jacoban (jaco set);
%  nob =  1: fixed-point iteration.
%  lib =  functions from HBVMlib1.
%  2020-02-13.
%__________________________________________________________________________
corH   = k>s;      % flag for RIE Hamiltonian corrections %
m      = length(y0);
y      = zeros(ns+1,m);
t      = (stride*h)*(0:ns)';
y(1,:) = reshape( y0, 1, m ); 
ierr   = 0;
itot   = 0;  
conta  = 0;
steps  = 0;
gs     = -1;
if ns<=0, return, end  % no output points %  
[c,beta,PO,Xs,Is] = lib.notree(k,s);        % HBVM(k,s) at Gaussian nodes %
beta = h*beta(:).'; c = h*c; Is = h*Is;
if nob==1 
   CC = []; teta = []; 
else
   I = eye(m); gamma = min(abs(eig(Xs))); CC = inv(Xs)*gamma;
end   
if k>s, H0 = feval( fun, y(1,:) ); end   % for RIE Hamiltonian correction %
ti    = t(1);
yi    = y(1,:);   
ii    = 1;
local = 0;
for i = 1:n
   local = local+1; 
   [LL,J0] = lib.computeL( nob, fun, ti, yi, h, jaco );
   if nob<=0, teta = inv( I -(h*gamma)*J0 ); end
   [yi,flag,it,tol,gsi] = ...                       % nonlinear iteration %
                   lib.itera( c, beta, PO, Is, fun, ti, yi, LL, teta, CC ); 
   if flag~=0
      warning('at step %i',i), ierr = ierr+1; 
   end
   if flag<0 
      warning('NaN: Exit'), y = y(1:ii,:); t = t(1:ii); ierr = -1; break
   elseif (k>s) && corH
      %ti = ti+h; [yi,conta] = lib.correggiH( fun, ti, yi, H0, tol, conta );
      ti = ti+h;
   else
      ti = ti+h;
   end
   if local==stride
      ii = ii+1; y(ii,:) = yi; local = 0;
   end   
   itot = itot +it;
   gs   = max(gs,gsi);
end 
steps = i;
return
 
function [y,t,itot,ierr,conta,steps,gs] = ...
                  bhbvm2( fun, y0, k, s, h, n, nob, stride, jaco, ns, lib )
%__________________________________________________________________________
%
%  Fixed stepsize implementation of the  HBVM(k,s) method for second order
%  Hamiltonian problems in the form q'' = M *nabla U(q),  with  M = M^T.
%  nob =  0: blended iteration, numerical Jacobian;
%  nob = -1: blended iteration, analytical Jacoban (jaco set);
%  nob =  1: fixed-point iteration.
%  lib =  functions from HBVMlib1.
% 2020-02-13.
%__________________________________________________________________________
corH   = k>s;      % flag for RIE Hamiltonian corrections %
m      = length(y0);
m2     = m/2;
iq     = 1:m2;     % indexes of the q-entries %
ip     = m2+1:m;   % indexes of the p-entries %
y      = zeros(ns+1,m);
t      = (stride*h)*(0:ns)';
y(1,:) = reshape(y0,1,m); 
ierr   = 0;
itot   = 0;  
conta  = 0;
steps  = 0;
gs     = -1;
if ns<=0, return, end  % no output points %  
[c,beta,PO,Xs,Is] = lib.notree(k,s);        % HBVM(k,s) at Gaussian nodes %
beta   = h*beta(:).'; c = h*c; Is = (h*h)*Is*Xs; w = (h*h)*Xs(1,:)*PO;
if nob==1 
   CC = []; teta = [];
else
   I  = eye(m2); gamma = min(abs(eig(Xs))); CC = inv(Xs)*gamma;
   CC = CC*CC; hg2 = (h*gamma)^2;    % blending for second-order problems %
end   
if k>s, H0 = feval( fun, y(1,:) ); end   % for RIE Hamiltonian correction %
M     = feval( fun );
if max(size(M))>1, flagM = 0; else, flagM = M==1; end
ti    = t(1);
yi    = y(1,:);   
ii    = 1;
local = 0;
if nob<=0, hg2M = hg2*M; end
for i = 1:n
   local = local+1; 
   [Lip,J0] = lib.computeL( nob, fun, ti, yi(iq), h, jaco );
   if nob==0, J0 = ( J0+J0.' )/2; end  % symmetrize the numerical Hessian %
   if nob<=0, teta = inv( I -hg2M*J0 ); end   
   [yi,flag,it,tol,gsi] = lib.iteras( ...           % nonlinear iteration %
             flagM, M, h, w, c, beta, PO, Is, fun, ti, yi, Lip, teta, CC ); 
   if flag~=0
      warning('at step %i',i), ierr = ierr+1; 
   end
   if flag<0 
      warning('NaN: Exit'), y = y(1:ii,:); t = t(1:ii); ierr = -1; break
   elseif (k>s) && corH
      ti = ti+h; 
      [yi,conta] = lib.correggiHs( fun, M, ti, yi, H0, tol, conta );   
   else
      ti = ti+h; 
   end
   if local==stride
      ii = ii+1; y(ii,:) = yi; local = 0;
   end   
   itot = itot +it;
   gs   = max(gs,gsi);
end
steps = i;
return
 
function [y,flag,iter,tol,gsi] = ...
                  itera( hc, hbeta, PO, Is, fun, t0, y0, L, teta, CC, now )
%
%  Blended or functional iteration for solving the discrete problem.
%  INPUT:
%      - hc      : stepsize h times the abscissae (column) vector; 
%      - hbeta   : stepsize h times the (row) vector of the weights;
%      - PO      : calP_s^T * Omega;
%      - Is      : h * calI_s;
%      - fun     : function evaluating the rhs of the equation (see above);
%      - (t0,y0) : current point;
%      - L       : local Lipschitz constant;
%      - teta    : (I-h*rho_s J0)^(-T), for the blended iteration;
%      - CC      : rho_s * X_s^(-1), for the blended iteration;
%      - now     : flag, set if no warning are to be issued.
%  OUTPUT:
%      - y       : new approximation;
%      - flag    : 0 if OK; 1, if tol not reached; -1 if NaNs occurred;
%      - iter    : number of nonlinear iterations;
%      - tol     : worst exit tolerance;
%      - gsi     : norm of the last Fourier coefficient.
%
%  2020.02.13.
%
if nargin<11, now = 0; end % no warnings are issued %
nob    = isempty(CC);
m      = length(y0);
tol0   = 8*eps*sqrt(m);    % 8*eps*m
tol    = 2*sqrt(L)*tol0;   % L*tol0
[k,s]  = size( Is );
fatt   = 10; if nob, fatt = 2*fatt; end
itmax  = fatt*s*m;
t      = t0+hc;
gamma  = zeros(s,m);
ey0    = ones(k,1)*y0;  % y0 is a row vector %
y      = ey0;
err    = inf;
f      = feval( fun, t, y );
scal   = diag(1./(1+abs(y0)));
for iter = 1:itmax
   if nob     % fixed-point iteration %
      gamma = PO*f;
   else       % blended iteration %
      eta   = gamma -PO*f;
      eta1  = CC*eta;  
      gamma = gamma -(eta1 +(eta-eta1)*teta)*teta;
   end   
   yold  = y;         
   y     = ey0 +Is*gamma;
   f     = feval( fun, t, y );
   err0  = err; 
   err   = norm((y-yold)*scal);  
   if ( (err<=tol) && (err>=err0) )||( err<=tol0 ), break
   elseif isnan(err), break, end  
end
gsi = norm(gamma(s,:)); 
y   = y0 +hbeta*f;            % new approximation %
if isfinite(err)
   flag = (iter>=itmax) && (err>tol);   
   if flag && (now==0)
      warning('Iteration convergence error = %g',err); 
   end
else
   if now==0, warning('NaNs occurred'); end, flag = -1; 
end   
return

function [y,flag,iter,tol,gsi] = iteras( flagM, M, ...
        h, h2e1XsPO, hc, hbeta, PO, h2IsXs, fun, t0, y0, L, teta, CC, now )
%
%  Blended or functional iteration for solving the discrete problem for
%  second order problems in the form q'' = M\nabla U(q), M = M^T.
%  INPUT:
%      - flagM   : flag set to 1 if M==1;
%      - M       : symmetric matrix at the rhs of the equation or scalar
%                  if M = alfa*I (see above);
%      - h       : stepsize;
%      - h2e1XsPO: h^2 * e_1^T * X_s * P_s^T * Omega;
%      - hc      : stepsize h times the abscissae (column) vector; 
%      - hbeta   : stepsize h times the (row) vector of the weights;
%      - PO      : calP_s^T * Omega;
%      - h2IsXs  : h^2 * calI_s * X_s;
%      - fun     : function evaluating the rhs of the equation (see above);
%      - (t0,y0) : current point. y0 = [q0 p0];
%      - L       : local Lipschitz constant;
%      - teta    : (I-h^2*rho_s^2*nabla^2U(q0))^(-1) for the blended 
%                              iteration for special second-order problems;
%      - CC      :  rho_s^2*X_s^(-2), for the blended iteration;
%      - now     : flag, set if no warning are to be issued.
%  OUTPUT:
%      - y       : new approximation. y = [q p];
%      - flag    : 0 if OK; 1, if tol not reached; -1 if NaNs occurred;
%      - iter    : number of nonlinear iterations;
%      - tol     : worst exit tolerance;
%      - gsi     : norm of the last Fourier coefficient.
%
%  2020.02.13.
%
if nargin<14, now = 0; end  % no warnings are issued %
if isempty(CC), nob = 1; else nob = 0; end  %% flag for fixed-point %%
m     = length(y0);
m2    = m/2;
tol0  = 8*eps*sqrt(m);    % 8*eps*m
tol   = 2*sqrt(L)*tol0;   % L*tol0  
[k,s] = size( h2IsXs );
fatt  = 10*( 1 +nob );
itmax = fatt*s*m;
t     = t0+hc;
gamma = zeros(s,m2);
q0    = y0(1:m2);
p0    = y0(m2+1:m);
q0p0  = ones(k,1)*q0 +hc(:)*(p0*M);
Q     = q0p0;
err   = inf;
f     = feval( fun, t, Q ); % *M is inside fun %
scal  = 1./(1+abs(q0p0));
for iter = 1:itmax
   if nob     % fixed-point iteration %
      gamma = PO*f;
   else       % blended iteration %
      eta   = gamma -PO*f;   eta1  = CC*eta;  % equivalent forms %
      gamma = gamma -(eta1 +(eta-eta1)*teta)*teta;
   end   
   Qold  = Q;         
   if flagM
       Q = q0p0 +h2IsXs*gamma;
   else    
       Q = q0p0 +h2IsXs*(gamma*M);
   end    
   f     = feval( fun, t, Q );   
   err0  = err; 
   err   = norm((Q-Qold).*scal);       
   if ( (err<=tol) && (err>=err0) )||( err<=tol0 ), break
   elseif isnan(err), break, end  
end
gsi = norm(gamma(s,:));
y   = [(q0+h*p0*M+h2e1XsPO*f*M) (p0+hbeta*f)];  % new approximation %
if isfinite(err)
   flag = (iter>=itmax) && (err>tol);   
   if flag && (now==0)
      warning('Iteration convergence error = %g',err); 
   end
else
   if now==0, warning('NaNs occurred'); end, flag = -1; 
end   
return

function [y,conta] = correggiH( fun, t, y, H0, tol, conta )
%
%  Hamiltonian correction, if of the order of 
%  Roundoff-Iteration Errors (RIEs).
%  INPUT:
%       - fun   : evaluates the ODE rhs (see above);
%       - t,y   : current point. y = [q p] or [p q];
%       - H0    : correct Hamiltonian value;
%       - tol   : exit tolerance of the nonlinear iteration
%                 which has computed y;
%       - conta : counter of points exceeding the RIE error.
%  OUTPUT:
%       - y     : updated point;
%       - conta : updated counter.
%
%  2020.02.13
%
if isempty(H0), return, end
H    = feval( fun, y );  % current Hamiltonian %
d    = feval( fun, t, y );      
nd   = d*d';
DH0  = tol*max(1,norm(y))*max(1,sqrt(nd));   
if abs(H-H0)>DH0  % Hamiltonian error too large, to be due to RIEs %
   conta = conta+1; return     
end    
m2   = length(y);
m    = m2/2;
d    = d(:).';
alfa = (H0-H)/nd;
d    = [-d(m+1:end) d(1:m)];    % the problem is in canonical form %
H1   = feval( fun, y+alfa*d );  % y = [q p] %
if abs(H0-H1) < abs(H-H0), y = y+alfa*d; end
return

function [y,conta] = correggiHs( fun, M, t, y, H0, tol, conta )
%
%  Hamiltonian  correction,  if  of  the  order  of 
%  Roundoff-Iteration Errors (RIEs), for problems
%  in the form  q'' = M*nabla^2 U(q),   M=M^T.
%  INPUT:
%       - fun   : evaluates the ODE rhs (see above);
%       - M     : symmetric matrix at the rhs of the equation or scalar
%                 if M = alfa*I (see above);
%       - t,y   : current point. y = [q p];
%       - H0    : correct Hamiltonian value;
%       - tol   : exit tolerance of the nonlinear iteration
%                 which has computed y;
%       - conta : counter of points exceeding the RIE error.
%  OUTPUT:
%       - y     : updated point;
%       - conta : updated counter.
%
%  2020.02.13
%
if isempty(H0), return, end
H    = feval( fun, y ); % current Hamiltonian %
m    = length(y);
iq   = 1:m/2;
ip   = iq+m/2;
d    = feval( fun, t, y(iq) );  
Mp   = y(ip)*M;
nd   = d*d' + Mp*Mp'; 
DH0  = tol*max(1,norm(y))*max(1,sqrt(nd));   
if abs(H-H0)>DH0        % Hamiltonian error too large to be due to RIEs %
   conta = conta+1; return      
end    
d    = d(:).';
alfa = (H0-H)/nd;
d    = [-d Mp];                 % the problem is in canonical form %
H1   = feval( fun, y+alfa*d );  % y = [q p] %
if abs(H0-H1) < abs(H-H0), y = y+alfa*d;  end
return

function [LL,J0] = computeL( nob, fun, ti, yi, h, jaco )
%
%  Approximate the Lipschitz constant and computes the
%  transpose of the Jacobian in case blended iteration 
%  is used. 
%  INPUT:
%       - nob   : flag for blended iteration (0) 
%                   or fixed-point iteration (1);
%       - fun   : evaluates the ODE rhs (see above);
%       - ti,yi : current point;
%       - h     : stepsize;
%       - jaco  : function evaluationg the Jacobian (see above)
%                 or void, if it has to be numerically approximated.
%  OUTPUT:
%       - LL    : approximated (local) Lipschitz constant;
%       - J0    : transpose Jacobian or void.
%
%  2015.02.14.
%
if nob==1  % fixed-point iteration %
   f0 = feval( fun, ti, yi ); 
   f1 = feval( fun, ti, yi+h*f0 );
   LL = 2*max( 5, norm( f1-f0 )/( h*max(1,norm(f0)) ) );
   J0 = [];                               % Jacobian not needed %
else       % blended iteration %
   if nob==0  
      J0 = numjacobian( fun, ti, yi ).';  % numerical Jacobian  %
   else       
      J0 = feval( jaco, ti, yi ).';       % analytical Jacobian %
   end   
   LL = max( 5, norm( J0, inf ) );
end
return

function Jf = numjacobian( fun, t0, y0 )
%  Auxiliary function for computeL.
%
%  NUMERICAL JACOBIAN  (slightly adapted from the code BiM/BiMD).
%  INPUT:
%       - fun   : evaluates the ODE rhs (see above);
%       - t,y   : current point.
%  OUTPUT:
%       - Jf    : approximated Jacobian.
%  2015.02.15
%
m  = length(y0);
Jf = zeros(m);
f0 = reshape( feval( fun, t0, y0 ), m, 1 );
dd = eps^(1/3);
for i = 1:m
  ysafe   = y0(i);
  delt    = sqrt( eps*max(dd,abs(ysafe)) );
  y0(i)   = ysafe +delt;
  f1      = feval( fun, t0, y0 );
  Jf(:,i) = ( f1(:) -f0 )/delt;
  y0(i)   = ysafe;
end
return

function [x,beta,PO,Xs,Is,Ps] = notree(k,s)
%
%  Gauss-Legendre weights, abscissae, and matrices for HBVM(k,s).
%
%  INPUT:
%       - k,s  : parameters for HBVM(k,s).
%  OUTPUT:
%       - x    : k Gauss-Legendre abscissae in [0,1];
%       - beta : corresponding quadrature weights;
%       - PO   : PO = P_s^T * Omega, with P_s(i,j) = P_{j-1}(x_i)
%                                    and Omega = diag(beta);
%       - Xs   : matrix X_s relating I_s and P_{s+1}: 
%                                    I_s == P_{s+1}[X_s; o^T];
%       - Is   : matrix I_s(i,j) = \int_0^{x_i} P_{j-1}(x)dx;
%       - Ps   : matrix P_s(i,j) = P_{j-1}(x_i);
%
%                i = 1...k, j = 1...s.
%
%  2020.02.13
%
[x,beta,xx] = lgwt01(k); 
[x,ind]     = sort(x); 
beta        = beta(ind); 
xx          = xx(ind);  
Ps          = zeros(k,s+1);  
for i = 0:s
    y = legendre(i,xx); Ps(:,i+1) = y(1,:).'*sqrt(2*i+1);
end   
j  = 1:s-1; xis = 0.5/sqrt(4*s^2-1);
Xs = (0.5)./sqrt(4*j.^2-1);
Xs = diag( Xs, -1 ); 
Xs = Xs-Xs.'; Xs(1,1) = .5;  
Is = Ps*[Xs; zeros(1,s-1) xis];
Ps = Ps(:,1:s);
PO = ( diag(beta)*Ps ).';  
return

function [x,w,y] = lgwt01(N)
% Auxiliary function for notree.
%  
% This script computes the Legendre-Gauss nodes (x) and weights (w)
% on the interval [0,1] with truncation order N.
% Also the nodes (y) on [-1,1] are computed.
% (Adapted from the file written by Greg von Winckel 2004.02.25).
%
% 2020.02.13
%
a  = 0; b = 1;  % normalization on the interval [a,b]==[0,1] %
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
    L(:,1) = 1;  L(:,2) = y;   %  within epsilon of old points          %
%    Lp(:,1) = 0; Lp(:,2) = 1;
    for k = 2:N1
        L(:,k+1) = ( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp = (N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    y0 = y;  y = y0-L(:,N2)./Lp;
end
x = (a*(1-y)+b*(1+y))/2;      % tranform from [-1,1] to [a,b] %    
w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2; % Compute the weights %
w = ( w + w(end:-1:1) )/2;              % and symmetrize them %
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%