function [V,x,iter] = solveValue(v0, x0, Pi);

global L M c alpha delta beta g CRIT lambda;

V = v0;
x = x0;
check = 1;
iter = 0;

while check > CRIT
    W = getW(V, x);
    nx = getX(W);
    nV = getNextV(W, nx, Pi);
    check = max( max(max(abs((nV-V)./(1+nV)))),  max(max(abs((nx-x)./(1+nx)))));
    V = lambda * nV + (1-lambda)*V;
    x = lambda * nx + (1-lambda)*x;
    iter = iter+1;
end