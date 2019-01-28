
function [nodes, Phi] = setupCheb(n, lo, hi)

%n = number of nodes
%lo = lower bound of the state space
%hi = upper bound of the state space

i = 1:n';
nodes = ((hi + lo)/2 + ((hi - lo)/2)*cos((n - i + .5)*pi/n))';

Phi = chebFuncs(nodes, n, lo, hi);
