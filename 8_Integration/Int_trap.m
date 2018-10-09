function [ returnVal ] = Int_trap(f, a, b, N )
%Int_mid: Calculate a definite integral by using the midpoint rule for N
%evenly spaces segments
%   f - Function to integrate
%   a - lower bound
%   b - upper bound
%   N - Number of segments

    assert(a < b);
    assert(floor(N) == N);
    assert(N > 0);

    %Setup Grid: 
    h = (b - a)/N;
    X = a:h:b;

    %Setup weights
    w = 2*ones(length(X),1);
    w(1) = 1;
    w(end) = 1;
    w=w*(h/2);

    %Evaluate function on grid:
    fv = f(X);

    %Sum...
    returnVal = fv*w;

end

