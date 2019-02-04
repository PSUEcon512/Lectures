
function [Pi, p] = getPi(p0);

global L M c alpha delta beta g;

options = optimset('TolX', 1e-12, 'TolFun', 1e-12);
[p, fval] = fsolve(@priceFOC, p0, options);

[g2, g1] = meshgrid(g);

ep1 = exp(g1 - p);
ep2 = exp(g2 - p');

Pi = M*(ep1./(1 + ep1 + ep2)).*(p - c);