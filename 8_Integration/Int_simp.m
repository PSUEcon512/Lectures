function [ returnVal ] = Int_simp(f, a, b, N )
%Int_mid: Calculate a definite integral by using the midpoint rule for N
%evenly spaces segments
%   f - Function to integrate
%   a - lower bound
%   b - upper bound
%   N - Number of segments

    assert(a < b);
    assert(floor(N) == N);
    assert(N > 0);
    
    %Make sure N is even:
    if (mod(N,2) == 1)
        N = N+1;
    end

    %Setup Grid: 
    h = (b - a)/N;
    X = a:h:b;

    %Setup weights
    w = 2*ones(length(X),1);
    evens = 2:2:N;
    w(evens) = 4;
    w(1) = 1;
    w(end) = 1;
    w=w*(h/3);

    %Evaluate function on grid:
    fv = f(X);

    %Sum...
    returnVal = fv*w;

end

