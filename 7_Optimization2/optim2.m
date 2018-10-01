%% Lecture 6 - Direction-Set and Quasi-Newton Optimization Methods
%
% Last week we examined optimization algorithms, ending with Newton's
% Method. We'll pick up there with quasi-Newton methods, which generally 
% replace the Hessian (or inverse Hessian) in the Newton Step with some approximiation. 
%

%% Stopping Criteria for Newton and Quasi-Newton Methods
%
% Unlike bracketing, or the area of the simplex for Nelder-Mead. We can't
% simply think about a shrinking interval for Newton methods. Instead, we
% use two stopping rules based on the information we've computed.
%
% First, has the sequence of guesses converged: 
%
% $$ ||x^{(k)} - x^{(k+1)}|| < \varepsilon(1 + ||x^{(k)}||)$
%
% This is a relative stopping rule, since it is scaled by the norm of
% $x^{(k)}$, the 1 simply deals with the case where $x^{(k)}$ is near the origin. 
% 
% Second, we want to know if we are at a local miniumum: 
%
% $$ ||\nabla f(x^{(k)}) || < \delta(1 + f(x^{(k)}))$$
%
% If this condition does not hold, then we have converged to a non-optimal
% point, and cannot claim the problem has been solved. Even if it does
% hold, we want to also check the first condition to make sure the problem
% doesn't just have a ``plateau'' near the optimum. 
%
% How to choose $\delta$ and $\varepsilon$?
%
% * Can't be lower (or same order as) machine epsilon, $eps$. 
% * Should not be tighter than the accuracy with which you are computing $\nabla
% f(x^{(k)})$ (finite differences?).
% * Often use values between $eps^{(1/4)}$ and $eps^{(1/2)}$. 

%% Overview: Generic Direction Set Methods
% Judd offers a ``generic
% algorithm'' which incorporates both a Hessian approximation and a line
% search: 
%
% *Initialization: Choose initial guess $x^{(0)}$ and stopping parameters
% $\delta$ and $\varepsilon > 0$. 
%
% # Compute search direction $s^k$. 
% # Find $\lambda_k$ that solves $\min_{\lambda} f(x^{(k)} + \lambda s^k)$. 
% # Assign $x^{(k+1)} \rightarrow x^{(k)} + \lambda_k s^k$. 
% # If $||x^{(k)} - x^{(k+1)}|| < \epsilon(1 + ||x^k||)$ go to 5, else, go
% to 1. 
% # If $||\nabla f(x^{(k)}) || < \delta(1 + f(x^{(k)}))$, Return convergence to local optimum,
% else, return convergence to non-optimal point. 
%
% The two ingredients are picking a direction, and a step length. First,
% let's consider how to pick a direction, keeping step length fixed. 
%
%% Steepest Descent
%
% The simplist quasi-Newton method simply replaces the Hessian with the
% identity matrix. This means the direction is entirely determined by the
% gradient. 
%
% $$ x^{(k+1)} \leftarrow x^{(k)} - \nabla f(x^{(k)}) $$
%
% Advantages and disadvantages: 
%
% * Always moves towards an optimum. 
% * Ignores curvature information
% * Linear convergence only. 
%
% In practice, this is not likely to work unless we dampen the step length.
% either by adding a line search or simply shortening the step size using a
% constant factor $a$. Here is an implemention edited from our Newton code last
% week:
%
% Recall we want to minimize:
%
% $$f(x)=(x_1-x_2)^4+2x_1^2+x_2^2-x_1+2x_2$$
%
% So we use:
n=0;            %initialize iteration counter 
focerr=1;          %initialize error 
x = [1; 1];
%x=[.05;-.5];        %set starting value
%
a = 0.09; %dampen value, if we set this anywhere near 1, things will blow up.

%Computation loop 
while focerr>1e-5&n<100 
    %Compute Gradient
    gradf=[4*(x(1)-x(2))^3+4*x(1)-1;...
          -4*(x(1)-x(2))^3+2*x(2)+2];                            
    
    %Compute Hessian (only for Newton)
    %Hf=     [12*(x(1)-x(2))^2+4,  -12*(x(1)-x(2))^2;...           
    %        -12*(x(1)-x(2))^2,    12*(x(1)-x(2))^2+2]; 
    
    %Perform Iteration
      %Newton: 
      %y=x-Hf\gradf 
      %Steepest Descent:
      y = x - a*gradf;
    
    x=y;                                                          
    n=n+1; 
    
    %Calculate FOC Error
    focerr= norm(gradf);                            
                                                          
end 
n,x,focerr,        %display end values

%%
% Recall that Newton's method solved this problem in only 8 iterations.
% Also, check what happens when we raise $a$ closer to 1. 

%% Davidson-Fletcher-Powell (DFP) and Broyden-Fletcher-Goldfarb-Shano (BFGS)
%
% Steepest descents biggest disadvantage is that it completely ignores
% curvature. DFP and BFGS are more akin to our nonlinear equations
% quasi-Newton methods in that they use consecutive iterates to estimate
% curvature. 
%
% Both satisfy the *quasi-Newton* Condition. Whereas a Newton step is approximately,  
%
% $$ x^{(k+1)} - x^{(k)} = d^{(k)} \approx  H^{-1}(x^{(k)}) [ \nabla f(x^{(k)} + d^{(k)} ) - \nabla f(x^{(k)})] $$
%
% Where the approximation comes from a finite difference intuition. The
% quasi-Newton condition selects a hessian approxmation to satisfies this exactly: 
%
% $$ d^{(k)} =  B^{(k)} [ \nabla f(x^{(k)} + d^{(k)} ) - \nabla f(x^{(k)})] $$
%
% DFP and BFGS are simply two different updating schemes that satisfy this
% condtion. They are both generalizations of the secant method.
%
% DFP: 
%
% $$ B^{(k+1)} \leftarrow B^{(k)} + \frac{dd^T}{d^T u} - \frac{B^{(k)} u u^T B}{u^T B^{(k)} u} $$
%
% where $u = \nabla f(x^{(k+1)}) - \nabla f(x^{(k)})$. 
%
% BFGS: 
%
% $$ B^{(k+1)} \leftarrow B^{(k)} + \frac{1}{d^T u} \left( wd^T + dw^T - \frac{w^Tu}{d^T u} dd^T \right) $$
%
% where $w = d - B^{(k)}u$. 
%
% I'm going to spare you the <https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm algebra>, suffice it to say that both extend the intuition of Broyden's method to Hessian approximation and BFGS solves
% the dual problem that DFP solves. In practice BFGS is the market leader. 
%
% Miranda and Fackler implement Steepest *Ascent*, DFP and BFGS updating in
% their maximization function |qnewton|.
%
% The MATLAB native function |fminunc| implements BFGS updating under the
% option |'Algorithm' = 'quasi-newton'| (which is also the default). 
%
% How does it work on our toy problem. First, let's writhe our function out
% explicilty:
%
% <include> qfunc.m </include>
%%
% Now we can call |fminunc|:
options = optimoptions('fminunc','Algorithm','quasi-newton',...
          'SpecifyObjectiveGradient',true, 'Display','iter');
[xstar, fstar] = fminunc('qfunc', [1, 1], options);
disp(xstar);

%%
% Can also use numerical differentiation, we get more function evaluations,
% but the same number of iterations:
options = optimoptions('fminunc','Algorithm','quasi-newton',...
          'SpecifyObjectiveGradient',false, 'Display','iter');
fminunc('qfunc', [1, 1], options);

%%
% That's everything I have for determining the direction, now let's see
% what the "step-size" column is all about...

%% Line Search
%
% Once we have settled on a direction, we may want to try step sizes other
% than 1 to optimize the objective as much as posible in the chosen
% direction. This is simply a 1 dimensional search problem. If the search
% direction is $s_k$ then we want to solve: 
%
% $$\min_{\lambda} f(x^{(k)} + \lambda s^k)$$
%
% We could use golden-search for this, but we don't typically because it
% spends too much time optimizing intermediate problems, only need a very
% approximate solution. 
%
% Hence, it is useful to use backtracking based on the *Armijo-Goldstein
% condition*, for an initial step size $\alpha_0$ and control parameters
% $\tau \in (0,1)$ and $c \in (0,1)$:
%
% # Set $t = -c s_k' \nabla f(x)$ and iteration counter $j = 0$
% # While $f(x) - f(x + \alpha_j s^k) > \alpha_j t$, increment $j$ and set
% $\alpha_j = \tau \alpha_{j-1}$.
% # Return $x + \alpha_j s^k$ as next iterate. 
%
% In words, backtrack until the improvement in $f$ is at least $c$ of that
% predicted by taking the linear Taylor approximation around $x$. 
%
% There are many alternative line search methods, (see Miranda and Fackler
% for several), however this captures the basic intuition. 
%
% While |fminunc| has limited flexibility in line search options, the
% Miranda and Fackler package |qnewton| is easier to customize, but I've
% found the MATLAB functions to be more robust. 