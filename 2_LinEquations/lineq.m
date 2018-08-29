%% Class 2: Linear Equations
% This week, we'll briefly review solviing linear equations. In fact,
% solving linear equations will be at the heart of a lot of what we do. For
% example Newtons method for solving a non-linear system can be thought as
% an iterative procedure that successively solves a linear approximate
% problem to the non-linear problem of interest. 
%
% Almost always, we want to solve non-linear systems using package code. So
% this class is mostly about making sure we have some idea what that code
% does.
%

%% Linear Equations
% We'll devote today to solving problems of the form
%
% $$ Ax = b $$
%
% This is easily the most common problem you will need to solve in
% empirical work, for a quick example OLS is an example of this problem
% where, $A = X'X$ and b = $X'y$. 
%
% To keep things simple, let's keep our matrices much smaller than we
% normally would: 
A =  [28    69    44    19
     5    32    38    49
    10    95    77    45
    82     3    80    65 ];

b = [100 ; 250 ; 300; 400];
%%
% You might think the way way to use a matrix inversion, after all we'd analytically
% solve this by writing $x = A^{-1}b$, and in fact this
% will work

x = inv(A)*b
%%
% However, this first computes the inverse (although maybe modern MATLAB
% catches this for you, that is still what you are asking for). In fact,
% computing the inverse takes many more computational steps than just
% solvng for $x$. 
%
% So what we really want is: 
x = A\b
%%
% For a little intuition, using backslash tells Matlab exactly what you
% want, rather than asking it to compute an intermediate object |inv(A)|
% which isn't what you are ultimately interested. 
%
% We don't need to worry much about how the computer solves linear
% equations because there are VERY good libraries that do it well. Most of
% the time we just call them. Today, we'll look a little more closely at
% what is going on. 

%% Direct Methods
%
% So what is MATLAB actually doing? As long as matrices are "small" (which
% is actually quite large these days). We can compute them directly using
% a standard series of steps.  This is great because it 
%
% * Is guarenteed to get us an answer
% * Completes in finite time.
%
% The downside is if the matrices are very large, finite time can still
% take longer than we want, which is why we will also talk about iterative
% methods. 

%%% L-U Factorization
%
% The key insight to LU factorization is solving a triangular system of
% equations is straightforward. 
%
% # Start with equation involving single variable. Solve it with scalar
% algebra. 
% # Now move on to equation with original variable plus additional, it is
% now an equation of a single variable, given solution in 1. Solve it. 
% # Continue until entire system is solved
%
% Because things need fancy names, lets call this _forward substiution_,
%
% To see how this works, lets use a triangular matrix.
T = tril(A) 
%%
% And solve
%
% $$ Tz = b$$ 
%
%using _forward substiution_ 
z = zeros(size(b))
z(1) = b(1)/T(1,1)
for i = 2:length(b)
    z(i) = (b(i) - T(i, 1:i-1)*z(1:i-1)) / T(i,i)
end
%%
%
% Produces the same result as: 
z = T\b
%
%%
% Great, but we want to be able to solve $Ax = b$ where $A$ is not
% triangular. If we can decompose $A$ into triangular matricies such that $A = LU$. We'd be
% golden. As it happens...

[L, U] = lu(A)

%% 
% Now we have $Ax = LUx = L(Ux)=b$, so we can
% 
% # Solve $Ly = b$ using _forward substitution_
% # Solve $Ux = y$ using _backward substitution_ (or just reorder rows and
% cols)

y = L\b
x = U\y
%%
% which is the same as: 
A\b
%%
% All that is left is to figure out how to do the decomposition, which is
% just an application of Gaussian Elimination. 

%%% Gaussian Elimination
%
% Remember elementary row operations on matrices from high school linear
% algebra? Neither do I. But all LU decomposition does is apply elementary
% row operations to convert $IA=A$ to $LU=A$
mL = eye(size(A));
mU = A;

%%
% I first want to get rid of the 5 in |mU(2,1)|, if I multiply the first
% row by 5/28 and add that to the second row, I'll eliminate the 5. Of
% course I'll need to account for this action in the $mL$ matrix.
%
mL(2,1) = mU(2,1)/mU(1,1)
mU(2,:) =  mU(2,:) - mU(1,:)*mL(2,1) 
mL*mU
%%
%
% Now we just need to loop over all the values we need to eliminate (let's
% start from scratch so the loop is nice): 

mL = eye(size(A));
mU = A;
nEq = length(b);
for c = 1:nEq-1
    for r = c+1:nEq
        mL(r,c) =  mU(r,c)/mU(c,c);
        mU(r,:) =  mU(r,:) - mU(c,:)*mL(r,c);
        %[mU; mL]; %Output Suprressed
    end
end
disp(mL);
disp(mU);
%%
% And now we know how to solve for x: 
y = mL\b;
mx = mU\y

x = A\b

%% Pivoting
%
% Notice, however that that the $L$ and $U$ I computed and those that
% |lu(A)| computed were different. Why is that? 
%
% * So far, we've been assuming
%   that a computer stores real numbers, in fact, it stores <https://en.wikipedia.org/wiki/Floating-point_arithmetic floating point
%   numbers> that consist of an exponent and significand (or mantissa). 
% * Since
%   the mantissa is of finite length, at some point we will need to round,
%   which can cause rounding errors. 
% * These isues can comound on themselves when we want to manipulate numbers of very different
%   magnitude. 
% * The MATLAB |lu(A)| (and other built in operations) include extra
%   _pivoting_
%   code to try to avoid this problem. Notice how its
%   decomposition has similar magnitude entries.
% * More or less, pivoting amounts to interchanging rows and columns to
% make the diagonal elements as large as possible, (since they are in the
% denominator). 

%%% Aside: Floating Point and its Discontents
% 
% Here is a quick example of a floating point pitfall:
for ind = 0:.1:1;
    if ind ~= fix(10*ind)/10
        disp(ind)
        disp(ind - fix(10*ind)/10)
    end
end
%%
% * Why does anything print out?
% * One moral of this story: *do not use logical equalities and floating point.* 
% * A second moral: Be careful when working with variables of very
% different magnitude, you may want to just rescale them so their magnitudes are similar. 
%
%% Conditioning, good and ill 
%
% * While pivoting can deal with rounding error issues in "nice" matrices,
% some matrices are not "nice" meaning they are inherently difficult, or
% impossible, to solve. 
%
% * The fact that such matrices exist is obvious, you can't solve a $Ax =
% b$ if $A$ is singular. A computer will have trouble if $A$ is "nearly"
% singular, but how to be more precise about what "nearly" means?
%
% * Following Judd, suppose we perturb $Ax = b$ to $A(x+e) = b + r$ so we
% want to know how big $e$ is for a small perturbation $r$. 
%
% * Linearity implies $Ae = r$ and in turn: 
%  $||A||\cdot||e|| \geq ||r||$, $e = A^{-1}r$,
%  and $||e|| \geq ||A^{-1}|| \cdot ||r||$. For any valid matrix norm.
%  Let's assume the $\ell_2$ norm for today.
%
% * Then let's define the elasiticty of error with respect to a
%  perturbation of b, $\frac{||e||}{||x||} \div \frac{||r||}{||b||}.$
%
% * It is interesting because it tells us how rounding errors in $b$ "grow"
%  through the solution of $Ax = b$. Bigger elasticity is going to cause more numerical problems. 
%
% * Now we can bound to percentage error: 
%
% $$ \frac{||r||}{||A||||A^{-1}|| ||b||} \leq \frac{||e||}{||x||} \leq
% \frac{||A||||A^{-1}|| ||r||}{ ||b||} $$
%
% And define the condition number as |cond(A)| $=||A||||A^{-1}||$ it
% gives an least upper bound on the error that is tight for some $b$. 
% It also gives us an idea of "almost singular": 

Z = [2 0 0
     0 1 0
     0 0 1e-8];
cond(Z)
norm(Z)
norm(inv(Z))
%%
% * The rule of thumb is that you lose one digit of precision for every
%  power of 10 in the condition number. 
%
% * MATLAB double precision floating
%  points are 64 bits: 1 sign, 11 exponent, 52 significant. This leads to
%  16 (base 10) digits of precision, $2^{-53} \approx 1.11 \times 10^{-16}$. 
%  MATLAB uses this rule of thumb to generate warning messages:

inv([2 0 0
     0 1 0
     0 0 1e-15])
%%
% Works fine, but, 
 
 inv([3 0 0
     0 1 0
     0 0 1e-16])
%%
% Gives us a warning. That said, we are losing a lot of precision in both cases.  
 
%% Iterative Methods
%
% So far, we've talked about direct methods of solving linear equations,
% and that is what we will use 99.9 percent of the time. We'll talk about
% iterative methods for two reasons
%
% # You want to work with **Really** big data, where the direct methods will be very time
%   consuming. If you find yourself here, you probably have data from Amazon,
%   Facebook, or similar. 
% # More likely, we're introducing iterative techniques in a linear context
% because they will be useful in more advanced context later. 
%
% Iterating is just repeatedly applying an operator in hopes that it converes to a
% fixed point which represents the solution of your problem.  
%
% In general, let $Q$ be some easy to invert matrix, then to solve $Ax = b$
% we can use the fact that 
%%
% $$Qx = b + (Q-A)x$$
%%
% $$x = Q^{-1}b + (I - Q^{-1}A)x$$
%
% So for some initial guess, $x^{(0)}$ we can iterate: 
%%
% $$x^{(k+1)} = Q^{-1}b + (I - Q^{-1}A)x^{(k)}$$
%
% If the process converges such that $||x^(k+1) - x^{(k)}|| < \epsilon$
% we've solved $Ax = b$. 
%
%
% The simplist possible version of this is $Q = I$: 
x = [1; 1; 1; 1;];
maxit = 100;
tol = 1e-4;
A = A .* ((10*eye(4))+1); %I'm cheating here, but we'll talk about why later...
for it = 1:maxit
    dx = (b - A*x);
    x = x + dx;
    if norm(dx) < tol 
        fprintf('Converged: x(1) = %.4f, iter = %d\n', x(1), it);
        break
    end
    if mod(it, 10)==0
        fprintf('it: %d \t norm(dx) = %.4f\n', it, norm(dx)); 
    end
end
%%
% As you can see, that doesn't work so well. The reason is that A is a long
% way from diagonal dominant. Specifically we can check whether this system
% is a contraction by checking
norm(eye(4) - A)
%%
% and see it is clearly not since the result is greater than 1. 
%%% Gauss-Jacobi
%
% Gauss Jacobi iterations are probably the most intuitive iteration method
% you could expect to work. In essence, for a guess of $x^{(t)}$ it guesses $x^{t+1}$ by 
% solving for $x_i$ in the $i$th equation using fixing $x_{-i}$.  In short:
% 
% $$ x_i^{(k+1)} = x_i^{(k)} + \frac{1}{a_{ii}} \left(b - \sum_{j \neq i} a_{ij} x_j^{(k)} \right)$$
%
% This is pretty easy to code: 
x = [1; 1; 1; 1;];
maxit = 100;
tol = 1e-4;
d = diag(A);
for it = 1:maxit
    dx = (b - A*x)./d;
    x = x + dx;
    if norm(dx) < tol 
        fprintf('Converged: x(1) = %.4f, iter = %d\n', x(1), it);
        break
    end
    if mod(it, 1)==0
        fprintf('it: %d \t norm(dx) = %.4f\n', it, norm(dx)); 
    end
end
%%
% Just to make sure everything works
disp(x)
disp(A\b)
%%
% With a cleaner, more matrix-style notation, we get the same thing: 
x = [1; 1; 1; 1;];
Q = diag(diag(A));
for it = 1:maxit
    dx = Q\(b - A*x);
    x = x + dx;
    if norm(dx) < tol 
        fprintf('Converged: x(1) = %.4f, iter = %d\n', x(1), it);
        break
    end
    if mod(it, 1)==0
        fprintf('it: %d \t norm(dx) = %.4f\n', it, norm(dx)); 
    end
end
%%
% Where we get exactly the same answer as before. 
%%
% We can see that this should work because now the discount factor of our
% contraction is, 
norm(eye(4) - Q\A)
%%
% Of course, if it were above one, it still _might_ work for some starting values. This condition is sufficient but not necessary for convergence.  

%%% Gauss-Seidel
%
% Gauss-Seidel iteration goes a little bit further. Once you have computed
% $x_i^{(k+1)}$, why not use it immediately instead of waiting to the next
% round? 
%
% $$ x_i^{(k+1)} = x_i^{(k)} + \frac{1}{a_{ii}} \left(b - \sum_{j = 1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j = i+1}^{n} a_{ij} x_j^{(k)} \right)$$
%
% If you stare at this, (maybe for a while) you can see we are using
% $x^{(k)}$ with the upper
% triangular portion of $A$ and $x^{(k+1)}$ with the lower triangular
% portion (with diagonal). And hence it is equivalent to using 
Q = tril(A) 
%%
% in the "general" iteration method. To see this, we can write:
%
% $$Qx = b - Ux$$
%
% $$Qx = b - (Q - A)x$$
%
% Thus we can implement Gauss-Seidel as: 
x = [1; 1; 1; 1;];
Q = tril(A);
for it = 1:maxit
    dx = Q\(b - A*x);
    x = x + dx;
    if norm(dx) < tol 
        fprintf('Converged: x(1) = %.4f, iter = %d\n', x(1), it);
        break
    end
    if mod(it, 1)==0
        fprintf('it: %d \t norm(dx) = %.4f\n', it, norm(dx)); 
    end
end

%%
%
% Notice it converges faster, which is intuitive because it updates more
% quickly, the catch is it may be less stable. In our case we are lucky and
% Gauss-Seidel does produce a contraction
norm(eye(4) - Q\A)
%%
% However, keep in mind this doesn't have to be the case, consider our
% original A matrix, which worked just fine using direct methods,
A =  [28    69    44    19
     5    32    38    49
    10    95    77    45
    82     3    80    65 ];

x = ones(4,1);
Q = tril(A);
norm(eye(4) - Q\A)
%%
% Doesn't look good for our hero Gauss-Seidel, but maybe it'll work out. In fact it isn't good: 
success = 0;
for it = 1:maxit
    dx = Q\(b - A*x);
    x = x + dx;
    if norm(dx) < tol 
        success = 1;
        fprintf('Converged: x(1) = %.4f, iter = %d\n', x(1), it);
        break
    end
end
if ~(success)
    fprintf('Convergence Failed! norm(dx) = %.4f', norm(dx));
end

%%
% For this reason, if you can, you probably should direct methods for solving linear equations whenever possible. Of
% course, we'll have a use for these iterative methods later on in
% non-linear problems. That brings us to next week...