
function fs = chebEval(s, c, n, lo, hi)

%s = vector of scalars at which to evaluate the function
%c = parameters which define the approximation
%n = number of grid points/basis functions used
%lo,hi = endpoints of the interval over which funciton is approximated.

T = chebFuncs(s, n, lo, hi);

fs = T*c;