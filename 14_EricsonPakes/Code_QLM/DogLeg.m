function  [X, info, perf] = DogLeg(fun,par, x0, opts)
%DogLeg  Dog Leg method for nonlinear system of equations
%    f_i(x) = 0 , i=1,...,n
%  where  x  is a vector,  x = [x_1, ..., x_n] . 
%  In the discussion we also introduce the function
%    F(x) = .5 * sum(f_i(x)^2) .
%  The functions  f_i(x)  and the Jacobian matrix  J(x)  (with 
%  elements  J(i,j) = Df_i/Dx_j ) must be given by a MATLAB
%  function with declaration
%            function  [f, J] = fun(x, par)
%  par  may be dummy.
%  
%  Call:
%      [X, info {, perf}] = DogLeg(fun,par, x0, opts)
%
%  Input parameters
%  fun  :  String with the name of the function.
%  par  :  Parameters of the function.  May be empty.
%  x0   :  Starting guess for  x .
%  opts :  Vector with five elements:
%          opts(1) = Initial trust region radius.
%          opts(2:5) used in stopping criteria:
%              ||F'||inf <= opts(2)                     or 
%              ||dx||2 <= opts(3)*(opts(3) + ||x||2)    or
%              ||f||inf <= opts(4)                      or
%              no. of iteration steps exceeds  opts(5) .
%
%  Output parameters
%  X    :  If  perf  is present, then array, holding the iterates
%          columnwise.  Otherwise, computed solution vector.
%  info :  Performance information, vector with 6 elements:
%          info(1:4) = final values of 
%              [||f(x)||inf  ||F'||inf  ||dx||2  Delta] 
%          info(5) = no. of iteration steps
%          info(6) = 1 :  Stopped by small  ||f(x)||inf
%                    2 :  Stopped by small  ||F'(x)||inf
%                    3 :  Stopped by small x-step
%                    4 :  Stopped by  kmax
%                    5 :  Problems, indicated by printout.  
%  perf :  (optional). If present, then array, holding 
%            perf(1,:) = values of  ||f(x)||inf
%            perf(2,:) = values of  ||F'(x)||inf
%            perf(3,:) = Radius of trust region,  Delta
%            perf(4,:) = values of  beta

%  Hans Bruun Nielsen,  IMM, DTU.  99.06.10

   %  Check and initialize
   [x n f J] = check(fun,par,x0,opts);
   thrc = max(20,n)*eps;    % For checking consistency
   g = J'*f;   ng = norm(g,inf);   ng2 = norm(g);   nf = norm(f,inf);
   delta = opts(1);   kmax = opts(5); 
   F = (f'*f)/2;
   Trace = nargout > 2;
   if  Trace
         X = x*ones(1,kmax+1);
         perf = [nf; ng; delta; 0]*ones(1,kmax+1);
       end 
   k = 1;   nu = 2;   stop = 0;   nx = opts(3) + norm(x);   beta = 0;
   nh = 0;   % added 04.05.12

   while  ~stop
     %  Check stopping criteria
     if      nf <= opts(4),  stop = 1;
     elseif  ng <= opts(2),  stop = 2; 
     elseif  delta <= opts(3)*nx,  stop = 3; 
     else    %  Find step
       alpha = (ng2/norm(J*g))^2;   a = -alpha*g;   na = alpha*ng2;
       [Q R] = qr(J);   y = Q'*(-f);
       D = abs(diag(R));
       si = find(D <= thrc*max(D));  nsi = length(si);
       if  nsi    % Singular.  Check consistency
         if  norm(y(si)) > thrc*F
           stop = 5;
           disp('Singular, non-consistent Newton equations.')
         else  % Find minimum norm solution
           p = ones(1,n);   p(si) = zeros(1,nsi);   p = find(p);
           RR = R(p,p);   b0 = [RR\y(p); zeros(nsi,1)];
           N = [RR\R(p,si); -eye(nsi)];
           b = b0 - N*(N\b0);
         end
       else,  b = R\y;  end
       if  ~stop    %  Proceed with Dog Leg
         nb = norm(b);
         if      nb <= delta    % Newton step
           h = b;   beta = 1;   nh = nb;   dL = F;
         elseif  na >= delta    %  Steepest descent
           h = -(delta/ng2)*g;   beta = 0;   nh = delta;
           dL = delta*(ng2 - .5*delta/alpha);
         else    % 'True' dog leg
           c = b - a;   cf = [c'*[c  2*a]  na^2-delta^2];
           beta = max(roots(cf));
           h = a + beta*c;   nh = delta;
           dL = .5*alpha*(1-beta)^2*ng2^2 + beta*(2-beta)*F;
         end
         if  nh <= opts(3)*nx,  stop = 3; end
       end
     end
     if  ~stop    % Perform step
       xnew = x + h;   
       [fn,Jn] = feval(fun, xnew,par);   Fn = (fn'*fn)/2;
       dF = F - Fn;
       if  (dL > 0) & (dF > 0)
         x = xnew;   nx = opts(3) + norm(x);
         F = Fn;  J = Jn;  f = fn;  nf = norm(f,inf);
         g = J'*f;   ng = norm(g,inf);   ng2 = norm(g);
         delta = delta / max(1/3, (1 - (2*dF/dL - 1)^3));   nu = 2;
       else
         delta = delta / nu;  nu = 2*nu;
       end
       k = k + 1;
       if  Trace
             X(:,k) = x;
             perf(:,k) = [nf; ng; delta; 0];  perf(4,k-1) = beta;
           end
       if  k > kmax,  stop = 4; end 
     end
   end
   %  Set return values
   if  Trace
     X = X(:,1:k);   perf = perf(:,1:k);
   else,  X = x;  end
   info = [nf  ng  nh  delta  k-1  stop];

% ==========  auxiliary function  =================================

function  [x,n, f,J] = check(fun,par,x0,opts)
%  Check function call
   sx = size(x0);   n = max(sx);
   if  (min(sx) > 1)
       error('x0  should be a vector'), end
   x = x0(:);   [f J] = feval(fun,x,par);
   sf = size(f);   sJ = size(J);
   if  sf(2) ~= 1 | sf(1) ~= n
       tx = 'f  must be a column vector of the same length as  x';
       error(tx), end
   if  sJ(1) ~= sf(1)
       error('row numbers in  f  and  J  do not match'), end
   if  any(sJ ~= n),  error('J  must be an n*n matrix'), end
%  Thresholds
   if  length(opts) < 5
       error('opts  must have 5 elements'), end
   if  length(find(opts(1:5) <= 0))
       error('The elements in  opts  must be strictly positive'), end