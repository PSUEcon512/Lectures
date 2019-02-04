
clear;
global L M c alpha delta beta g CRIT lambda;

L = 18;
M = 5;
c = 5;
alpha = 3;
delta = 0.7;
beta = 0.925;
CRIT = 1e-10; %10^(-6);
lambda = 1 %.75; %For Dampening

g = zeros(L,1);
g(1:5) = 3.*[1:5]' - 4;
%for w = 6:L
%    g(w) = 12 + log(2 - exp(16 - 3*w));
%end
g(6:L) = 12 + log(2 - exp(repmat(16, L-5,1) - 3*[6:L]'));