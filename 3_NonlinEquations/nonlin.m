%% Class 3: Nonlinear Equations
%
% These notes will cover the basics of solving nonlinear equations. 
% Systems of nonlinear equations 
% show up in economics in a variety of ways, most obviously in the
% calculation of equilibrium. 
%
%% Nonlinear Equation Form 
%
% * Consider a function $f: R^n \rightarrow R^n$ our objective is to solve
% a rootfinding problem, that is, to deterimine $x$ such that $f(x) = 0$. 
%
% * We may also want to solve an
%   equivalent fixed point problem $g(x) = x$, simply by defining 
%   $f(x) = g(x) - x$.
%
%% Bisection
%
% * The simplist method of equation solving is bisection, which is applied
% only to single variable equations, although it can be used as part of
% multi-dimensional solution approaches. 
%
% * Bisection is based on the intermediate value theorem, so it requires:
%
% # $f : R \rightarrow R$ is continuous. 
% # There exists $a < b$ such that $f(a) < 0 < f(b)$. 
%
% * From here, the bisection algorithm simply bisects [a b] to get closer
% and closer to a zero: 
%
% # Given $[a, b]$, compute the midpoint $c \leftarrow (a + b) /2$ and $f(c)$.
% # if $f(c) < 0$, assign $a \leftarrow c$, else assign $b \leftarrow c$.
% # if $b - a < \epsilon$, stop; else go to 1. 
%
%
% Let's use it to solve something. First, define the function,
X = -2:.1:4;
Fx = 2 + exp(X) - 3.*(X.^2);
plot(X, Fx, X, zeros(size(X)));
%%
% It's continuous, and it has a zero (in fact, three) so bisection should be able to find a root. 
% Here's our bisection code: 
%
% <include> bisection.m </include>
%
%
%% 
% We still need to call it, which requires passing a function. In MATLAB,
% all variables are really pointers, so we can assign a function using
% anonomyous functions rather than writing a file for it: 
f = @(x) 2 + exp(x) - 3.*(x.^2);
f
f(0)
feval(f, 0)
%%
% Now we can call bisection: 
sol1 = bisection(f, -2, 4)

f(sol1)

%%
% Of course, we just found a solution, and my initial bracket actually
% contains three, a different bracket might find a different solution...
sol2 = bisection(f, -1, 4)

f(sol2)

%% 
% We can also adjust the stopping tolerance away from the defaults: 

sol3 = bisection(f, -1, 4, 1e-3, 1e-3)

f(sol3)

%% Iteration 
%
% As with linear equations, we can simply iterate nonlinear equations in
% hopes that they converge to a fixed point:
%
% $$x^{(k+1)} \leftarrow g(x^{(k)})$$
%
% This will work if we are sufficiently close to a fixed point $x^*$ where
% $||g'(x^*)|| < 1$.
%
%
g = @(x) [ x(1)^.5,  x(1).^.25 + x(2)^.75];

x = [.5 .5];

maxit = 1000;
for iter = 1:maxit
    nextX = g(x);
    if abs(x - nextX) < 1e-6
        break;
    end
    x = nextX;
    disp([iter x]);
end

format long
disp([iter x g(x)]);
format short
%%
% Of course, this is going to be far less robust than something like
% bisection, but it is easy, and if you know you have a contraction, it may
% be worth implementing before moving on to other methods. 
%% Newton's (aka Newton-Raphson) Method
%
% Newton's method is the real workhorse for solving nonlinear equations. It
% is an iterative scheme that follows the principle of _successive
% linearization_. That is, it (1) approximates the nonlinear problem with its
% linear taylor approximation, (3) solves that linear problem, and then (3) checks the solution
% of the original problem. If it is solved, we are done, if not, we
% iterate.
%
% Newton's method iterations are derived from a first-order Taylor approximation: 
%
% $$ f(x) \approx f(x^{(k)}) + f'(x^{(k)})(x - x^{(k)}) = 0 $$
%
% Solve for $x^{(k+1)}$ to get the iteration rule: 
%
% $$ x^{(k+1)} \leftarrow x^{(k)} - [f'(x^{(k)})]^{-1}f(x^{(k)}) $$
%
% For Newton's Method to work we need: 
%
% # $f$ is continuously differentiable.  
% # $f(x^{(0)})$ is sufficiently close to a root. 
% # $f'$ is invertible at the root.
%
% These will not always hold, and there is no theoretical way to determine
% sufficiently close. 
%
% Let's follow the Miranda and Fackler textbook and use Newton's method to
% solve a Cournot Duopoply. The primatives are: 
%
% * Demand for the good is CES: $P(q) = q^{-1/\eta}$, where $q = q_1 + q_2$. 
% * Each firm has quadratic costs: $C_i(q) = \frac{1}{2} c_i q_i^2$ for $i = 1,2$. 
%
% Then firm $i$ solves, 
%
% $$ \max_{q_i} P(q_i + q_{-i})q_i - C_i(q_i) $$
%
% We can find equilibrium of this model by simultaneously solving both
% firms' first order conditions, which will be nonlinear equations for $i = 1,2$: 
%
% $$ f_i(q) = (q_i + q_{-i})^{-1/\eta} - (1/\eta)(q_i + q_{-i})^{-1/(\eta-1)}q_i
% - c_iq_i = 0$$
%
% So we can write a function to calculate this function and its Jacobian: 
%
% <include> cournot.m </include>
%
%
% To solve this with Newton's method we just follow the iteration rule: 
q = [2; 2];
tol = 1e-8;
maxit = 100;
for iter =1:maxit
    [f, dF] = cournot(q);
    q = q - dF\f;
    if norm(f) < tol
        break
    end
end

fprintf('Iter: %d, q = (%f, %f), f(q) = (%f, %f)\n', iter, q(1), q(2), f(1), f(2));
%%
%
% The catch with Newton's method is that you have to compute the Jacobian,
% and things may not work very well if it is ill-conditionioned.
% This is simply because the Newton itertion solves a linear equation, and
% this can cause all the numerical issues we discussed last week.
% 
% Quasi-Newton methods, which we'll discuss next time, offer alternatives
% based on the successive linearization principle. 
%
% The MATLAB quasi-Newton method solver is fsolve:

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',true);
options
fsolve('cournot', [2; 2], options)
%%
% Seems to give us the same answer, that's good news. 
%
% By inspecting the steps, we can see one of the reasons that Newton's
% method is so attractive. Once we are in a neighborhood of the solution,
% Newton's method converges *quadratically* to the solution. 
%
% Of course, if you are not in the neighborhood, all bets are off. 

%% Numerical Differentiation
%
% While it is best to have an analytic derivative, it is common to
% approximate derivatives numerically. With most Newton-based solvers such
% as KNITRO and |fsolve| they will take their own numerical derivatives if
% you don't specify the gradient. In our example above, that also works
% just fine: 

options = optimoptions('fsolve','Display','final','SpecifyObjectiveGradient',false);
fsolve('cournot', [2; 2], options)

%% 
% What's going on under the hood? Exactly what you would expect. Recall the
% definition of a derivative:
%
% $$ f'(x) = \lim_{h \rightarrow 0} \frac{f(x+h) - f(x)}{h} $$
%
% So let's just approximate this with a 'small' $h$. How small is small?
% Back to Taylor expansions: 
%
% $$f(x + h) = f(x) + f'(x)h + O(h^2) $$
%
% $$ f'(x) = \frac{f(x+h) - f(x)}{h} + O(h)$$
%
% So the error is of order $h$. This is a one-sided finite difference
% equation. It is probably the most common because it requires only a
% single additional computation of $f$ for each direction. 
%
% However, its not hard to increase accuracy. Consider two second order
% expansions: 
%
% $$f(x + h) = f(x) + f'(x)h + f''(x)\frac{h^2}{2} + O(h^3) $$
%
% $$f(x - h) = f(x) - f'(x)h + f''(x)\frac{h^2}{2} + O(h^3) $$
%
% Subtracting one from the other and solve for $f'(x)$: 
%
% $$f'(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2)$$
% 
% We can take this further and gain even more accuracy at the cost of
% additional funciton evaluations. 
%
% This still leaves us the question of how to set $h$. Miranda and Fackler
% suggest the rule of thumb for one-sided derivatives. 
%
% $$h = \max(|x|, 1) \sqrt(\epsilon)$$
%
% Where $\epsilon$ is machine epsilon. The same rule of thumb for two-sided
% derivatives would lead to: 
%
% $$h = \max(|x|, 1) \sqrt[3](\epsilon)$$
%
% These amount to about $10^{-8}$ or to $10^{-6}$ respecitively.
%
% We can implement a one-sided numerical Jacobian as: 
%
% <include> myJac.m </include>
%
%
% Note that the Jacobian takes $N$ function evaluations to compute,
% therefore for large functions, computing a numerical derivative may
% produce a time sink if $f$ is expensive to compute.
% In this case, it may be worthwhile to code the analytical
% Jacobian.  This is a classic *programmer-time/compute time tradeoff*.  
%
% To see that it works: 
mJ = myJac('cournot', q)
[f, J] = cournot(q)
%%
%
%
% Now we can use it inside Newton's Method:
q = [2; 2];
tol = 1e-8;
maxit = 100;
for iter =1:maxit
    f = cournot(q);
    dF = myJac('cournot', q); 
    q = q - dF\f;
    if norm(f) < tol
        break
    end
end

fprintf('Iter: %d, q = (%f, %f), f(q) = (%f, %f)\n', iter, q(1), q(2), f(1), f(2));