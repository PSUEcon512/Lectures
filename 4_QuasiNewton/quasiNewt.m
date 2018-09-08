%% Class 4: Quasi-Newton Methods for solving Nonlinear Equations
%
% * Last week we introduced some methods for solving nonlinear equations, the
% workhorse method was Newton's method, which relied on iteratively
% linearlizing the nonlinear problem around an iterate. Leading to the
% iteration rule:
%
% $$ x^{(k+1)} \leftarrow x^{(k)} - [f'(x^{(k)})]^{-1}f(x^{(k)}) $$
%
% * One issue with Newton's method is that it required computing the Jacobian
% (matrix of derivatives) of the problem at every iteration, $[f'(x^{(k)})]^{-1}$. 
% While this can
% be done numerically it may be computatitionally intensive. 
%
% * Quasi-Newton methods are simply approaches to approximate the jacobian
% rather than computing it directly (either numerically or analytically). 

%% Secant Method
%
% Suppose $f : R \rightarrow R$ is univariate. 
%
% While computing the numerical derivative would only require a single
% function evaluation, we can save even that by just using the previous
% iterate: 
%
% $$f'(x^{(k)}) \approx \frac{ f(x^{(k)}) - f(x^{(k-1)}) }{ x^{(k)} -
% x^{(k-1)}}$$
%
% No the secant method iteration just replaces the derivative in Newtons
% method with this secant approximation: 
%
% $$x^{(k+1)} \leftarrow x^{(k)} - \frac{ x^{(k)} - x^{(k-1)} }{ f(x^{(k)})
% - f(x^{(k-1)}) }$$
%
% This illustrates the main concept of quasi-newton methods: use previously
% computed information to efficiently approximate the derivative
% information of the current iterate. 
%
% Formally, the secant method requires 2 initial guesses. Although usually
% we just compute the derivative for the first iteration. (e.g., use $x$ and $x+h$ as the initial guesses).  

%% 
% Let's recall our univariate function:
f = @(x) 2 + exp(x) - 3.*(x.^2);
X = -2:.1:4;
%Fx = 2 + exp(X) - 3.*(X.^2);
plot(X, f(X), X, zeros(size(X)))

%% 
% So to implement the secant method:

%Assign initial values
x = 0;
xOld = 1;
fOld = f(xOld);

% Secant iterations:
tol = 1e-8;
maxit = 100;
for iter =1:maxit
    fVal = f(x);
    fprintf('iter %d: x = %.8f, f(x) = %.8f\n', iter, x, fVal);
    if abs(fVal) < tol
        break
    else
        xNew = x - ( (x - xOld) / (fVal - fOld) )* fVal;
        xOld = x;
        x = xNew;
        fOld = fVal;       
    end
end


%% Broyden's Method
%
% Broyden's method generalizes the secant method to a multidimensional
% problem. However, how do we update the entire Jacobian with a single pair
% of function evaluations. 
%
% Effectively, for two values of $f:R^n \rightarrow R^n$ we need to solve 
%
% $$ f(x^{(k)}) - f(x^{(k-1)}) = A (x^{(k)} - x^{(k)}) $$
%
% But $A$ has $n^2$ values, and we only have $n$ equations. 
%
% We will supplement this with information from the previous approxmiation
% of the Jacobian, which we hope is "close" to the current Jacobian.
% Where close is determined according to the
% <https://en.wikipedia.org/wiki/Matrix_norm#Frobenius_norm Frobenius norm>.
%
% This leads to the iteration rule: 
%
% $$A^{(k+1)} \leftarrow A^{(k)} + [ f(x^{(k+1)}) - f(x^{(k)}) - A^{(k)}d^{(k)}] $$
% $$\frac{ d^{(k)^T} }{ d^{(k)^T} d^{(k)} } $$
%
% Where $d^{(k)}$ is the column vector $x^{(k+1)} - x^{(k)}$. Note that the
% term in brackets is $n \times 1$ while the fraction is $1 \times n$. 
%
% We could use $A^{(k+1)}$ directly, but then we would still need to solve
% a linear equation to get the next iterate. Instead, it turns out a matrix
% multiplication can directly deliver us an approximation of the inverse
% Jacobian. Then we have: 
%
% $$ B^{(k+1)} \leftarrow B^{(k)} + [ (d^{(k)} - u^{(k)}) d^{(k)^T}B^{(k)}]/(d^{(k)^T} u^{(k)})$$
%
% Working this out will take some algebra, but at least note that the
% denominator here is a scalar, so it is just a matrix multiplication in
% practice. 
%
% Let's implement it to solve our cournot model from last time. We'll need
% an initial guess at the Jacobian, so we'll use a numerical derivative:
q = [2; 3];
fVal = cournot(q)
iJac = inv(myJac('cournot', q))
%%
%
% Now for the Broyden iterations: 
maxit = 100; 
tol = 1e-6; 
for iter = 1:maxit
    fnorm = norm(fVal);
    fprintf('iter %d: q(1) = %f, q(2) = %f, norm(f(x)) = %.8f\n', iter, q(1), q(2), norm(fVal));
    if norm(fVal) < tol
        break
    end
    d = - (iJac * fVal);
    q = q+d;
    fOld = fVal;
    fVal = cournot(q);
    u = iJac*(fVal - fOld);
    iJac = iJac + ( (d - u) * (d'*iJac) )/ (d'*u);
end

%% 
% You will recall that using the derivative information directly took fewer
% iterations to converge, however each iteration was more computationally
% intensive since 
% 
% # Had to compute numerical Jacobean (or supply analyticl
% # Had to solve a linear equation as part of iteration. 

%% What's going wrong? 
%
% Inevitably, when you are trying to solve a nonlinear system, its not
% going to work (at least the first few attempts). Just remember... 
% *it's all your fault*.
%
% # Using a packaged solver can minimize coding errors in the solution
% algorithm itself. 
% # If you are coding your Jacobian, its a good idea to at least check your
% code against a numerical derivative. 
% # Check the coding of your function by computing it at some points where
% you can calculate the answer with paper and pencil. 
% # "Explore" (plot or grid search) your algorithm to attempt to find a
% good start point. 
% # Re-scale your function to avoid ill-conditioning. Try to keep
% "reasonable inputs" in the same order of magnitude. 
% # Be mindful of bounds, if you have a $log(x)$ in your equations, you
% defintely don't want to evaluate at $x = -2$, but the solver won't
% realize that unless you tell it. 
% # If you try to solve a system with no solution, the computer is not
% going to tell you this, it will just keep trying.
% # If you try to solve a system with kinks or discontinuities, your mileage
% may vary (to put it mildly).
% # Finally, Newton's method can always blame you for not being in the
% neighborhood of the solution. 
%
% Some last advice, sometimes you can transform your equations to make them closer to
% linear. In the extreme why solve: 
%
% $$ \exp(x) - 10 = 0 $$   
%
% When you can solve: 
%
% $$ x = log(10)$$
%
% However it may also be that a system like:
% 
% $$ x^{0.2} + y^{0.2} - 2 = 0 $$
%
% $$ x^{0.1} + y^{0.4} - 2 = 0 $$
%
% Can be more easily be solved after re-scaling to get closer to CRTS:
% 
% $$ (x^{0.2} + y^{0.2})^5 - 32 = 0 $$
%
% $$ (x^{0.1} + y^{0.4})^4 - 16 = 0 $$
%
% The broad lesson here is that the computer wants to solve the math
% problem you give it, if that problem needs to be manipulated slightly in
% ways that are not "economically intuitive" that is fine. In general, the
% closer to a _linear_ problem you have, the easier it will be to solve
% with Newton or quasi-Newton methods, since both follow the principle of
% successive linearization. 