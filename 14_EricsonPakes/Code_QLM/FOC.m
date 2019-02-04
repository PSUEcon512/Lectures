%
% Quality-ladder.
% FOC.
%

function [Del,Jac] = FOC(p,g);

% Globals.
global c;

% Demand.
u = exp(g-p);
d = u./(1+u);

% FOC.
Del = 1-(1-d).*(p-c);

% Jacobian.
Jac = -(1-d).*(d.*(p-c)+1);
