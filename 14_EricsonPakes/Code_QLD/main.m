clear;
tic;
setupParams;
[Pi, p] = getPi(repmat(c, L, L));
%load profit;

V0 = Pi./(1-beta);
x0 = zeros(L,L);
[Va, xa, ita] = solveValue(V0, x0, Pi);

elapsed = toc./60;
disp(sprintf('Elapsed time: %12.4f minutes', elapsed));

figure(1);
mesh(Va);
title('Value Function');

figure(2);
mesh(xa);
title('Policy Function');