%
% Quality ladder.
% Setup.
%

% Constants: Discounting.
beta = 0.925;				% discount factor.

% Transition probabilities: Quality states.
delta = 0.7;				% rate of depreciation.
alpha = 3;                  % effectiveness of investment.

% Constants: States.
L = 18;						% # states.

% Constants: Product market equilibrium.
w = [1:L]';
g = (3.*w-4).*(w<=5)+(12+log(2-exp(16-3.*w))).*(w>5);       % utility.
M = 5;					    % market size.
c = 5;                      % marginal cost.

% Constants: Program control.
tol = 1e-4;					% tolerance.
maxiter = 5e2;              % maximum number of iterations.