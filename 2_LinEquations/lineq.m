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

%%% Conditioing 

%% Iterative Methods

%%% When would we use them? 

%%% Gauss-Jacobi

%%% Gauss-Seidel