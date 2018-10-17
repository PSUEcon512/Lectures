function llf = llk(X,par,m)
% Calculate the (negative) log-likelihood function
% INPUTs:
%     X.a, X.b:  regressors (matrices)
%     X.y:       plant id and dependent variable (2-column matrix)
%     X.N:       number of plants (scalar)
%     par:       vecter of parameters
%     m:         number of points for Gaussian quadrature (1-by-2 matrix)
% OUTPUTs:
%     llf:       negative value of log-likelihood function

% Extract the coefficient parameters
p.a = par(1:19);                        % Coefficients of (a) equation
p.b = par(20:35);                       % Coefficients of (b) equation
p.w = [par(36) 0; par(37) par(38)];
p.w = p.w*p.w';                         % Omega matrix

% Get the points and associated weights for Gaussian quadrature
[x,w] = qnwnorm(m,[0 0],p.w);
k = length(w);

y = repmat(X.y(:,2),1,k);
m = size(y,1);
a = repmat(X.a*p.a,1,k)+repmat(x(:,1)',m,1);
b = repmat(X.b*p.b,1,k)+repmat(x(:,2)',m,1);
ind = repmat([0 0 0 1 1 1 1 1 1]',X.N,k);
phi = normcdf(a.*ind+b.*(1-ind));
phi = reshape(prod(reshape(phi.^y.*(1-phi).^(1-y),9,[])),X.N,[]);
llf = -sum(log(phi*w));
