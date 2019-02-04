%
% Qualiy ladder.
% Transition probabilities.
%

function p = transprob(w,x);

% Globals.
global L delta alpha;

% Transition probabilities.
p = zeros(1,L);
if w==1
    p(w+1) = ((1-delta).*alpha.*x(1))./(1+alpha.*x(1));
elseif w==L
    p(w-1) = delta./(1+alpha.*x(1));
else
    p(w+1) = ((1-delta).*alpha.*x(1))./(1+alpha.*x(1));
    p(w-1) = delta./(1+alpha.*x(1));
end
p(w) = 1-sum(p);
