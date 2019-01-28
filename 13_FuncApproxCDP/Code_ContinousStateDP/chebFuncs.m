
function T = chebFuncs(s, n, lo, hi)

%Convert the points to evaluate (s) into [-1,1]
z = (2*(s - lo)/(hi - lo)) - 1;

%Use recursive formula to generate the chebyshev polynomials. 
% Could use arccos formula, I do this just to make it clear these are
% polynomials...
T = zeros(length(s), n);
T(:,1) = ones(length(s), 1);
T(:,2) = z;
for j = 3:n
    %T(:,j) = 2*z.*T(:,j-1) - T(:,j-2);
    T(:,j) = 2*bsxfun(@times, z, T(:,j-1)) - T(:,j-2);
end