%% Lecture 8 - Numerical Integration
%
% For the next two weeks, we'll consider how to evaluate the definite integral of a
% real valued function over some interval $I \in R^d$:
%
% $$ \int_{x \in I} f(x) w(x) dx $$
%
% The weight function $w(x)$ will most commonly be a probability density function,
% in which case the integral denotes the expectation of $f(x)$. 
%
% There are plenty of places where such a problem will pop up:
%
% * Evaluating expectations of future states in a dynamic model. 
% * Integrating out unobserved heterogeneity in a likelihood estimation. 
% * Dealing with unobserved consumer types in a discrete choice or demand
% problem. 
% * ...
%
% All numerical integration formulae must break this problem into a finite
% number of calculations, so they all end up using approximations of the
% form: 
%
% $$ \int_{x \in I} f(x) w(x) dx \approx \sum_{i=0}^n w_i f(x_i)$$
%
% Where $w_i$ are the quadrature weights and $x_i$ are the quadrature
% nodes. The difference between these approaches is how the nodes and
% weights are chosen: 
%
% * Newton-Cotes methods approximate $f$ between nodes using polynomials. 
% * Gaussian quadrature choses nodes and weights which satisfy moment
% conditions.
% * Monte Carlo and quasi-Monte Carlo methods use "random" or
% equidistributed nodes and rely on asymptotic properties of their
% approximation. 

%% Newton-Cotes Formulae
%
% The simplest N-C formula is the *midpoint rule*, which simply uses the
% value of the function at the midpoint as its estimate of the function
% over the interval, it is formaly implied by: 
%
% $$ \int_a^b f(x) dx = (b - a) f\left(\frac{a+b}{2} \right) +
% \frac{(b-a)^3}{24}f''(\xi) $$
% 
% Where $\xi \in [a, b]$, and the sceond term represents the error in the
% midpoint approximation. Notice that the rule's error is increasing in the
% width of the interval. We can improve by making it a composit rule. 
%
% Suppose we divide $[a,b]$ into $n$ intervals of width $h = \frac{(b-a)}{n}$.
% Then we can write the rule: 
%
% $$ \int_a^b f(x) dx = h \sum_{j=1}^n f(x_j) + \frac{h^2(b-a)}{24}f''(\xi) $$
%
% Where $x_j = a + (j- \frac{1}{2})h, j = 1,2,\ldots,n$. Now the error is a
% function of $h^2$ so doubling the number of points will reduce the error
% by about 75 percent. 
%
% The midpont rule effectively approximates $f(\cdot)$ over an interval with a constant,
% aka, an order-0, polynomial. The *trapezoid rule* extends the logic to a
% line, an order-1 polynomial. 
%
% $$ \int_a^b f(x) dx = \frac{(b - a)}{2} [f(a) + f(b)] - 
% \frac{(b-a)^3}{12}f''(\xi) $$
%
% Now when we break the interval up, the endpoints will be treated differently than interior points, 
%
% $$ \int_a^b f(x) dx = \frac{h}{2} [f(x_0) + 2f(x_1) + \cdots + 2f(x_{n-1})
% + f(x_n)] - \frac{h^2(b-a)}{12}f''(\xi) $$
%
% where $x_i = a + ih$. Note that this fits into our weighted sum notation
% where $w_0 = w_n = h/2$ and $w_i = 2$ for $i = h,\ldots,n-1$. 
%
% This is worth coding up:
%
% <include> Int_trap.m </include>
% 
% And we can evaluate it: 
Int_trap(@(x) exp(x), 0, 1, 4)
%%
%Of course this integral is easy to compute analytically
%
% $$\int_0^1 \exp(x) dx = exp(1) - 1 $$
%
% So we can see the error as we add points: 
format long;
Int_trap(@(x) exp(x), 0, 1, 4) - (exp(1) - 1)
Int_trap(@(x) exp(x), 0, 1, 8) - (exp(1) - 1)
Int_trap(@(x) exp(x), 0, 1, 16) - (exp(1) - 1)
Int_trap(@(x) exp(x), 0, 1, 32) - (exp(1) - 1)


%%
% We can go further by approximating the function with piecewise quadratic
% funcitons instead of piecewise linear, which leads to *Simpson's Rule*. To fit a quadratic on the interval
% we need three points, so we bring back the midpont: 
%
% $$ \int_a^b f(x) dx = \frac{(b - a)}{6} [f(a) + 4f\left(\frac{a+b}{2} \right) + f(b)] - 
% \frac{(b-a)^5}{2880}f''''(\xi) $$
%
% Exactly why this is the formula for the area under a parabola may not be immediately obvious, but it
% is just some basic calculus. If you don't belive me this
% <https://www.youtube.com/watch?v=7MoRzPObRf0 youtube video> derives it. 
%
% We can then create a composite version of simpsons rule the only
% complication is that now we will use half size intervals to account for
% the midpoints. 
%
% $$ \int_a^b f(x) dx = \frac{h}{3} [f(x_0) + 4f(x_1) + 2f(x_2) + 4f(x_3) +  \cdots + 4f(x_{n-1})
% + f(x_n)] - \frac{h^4(b-a)}{180}f''''(\xi) $$
% 
% Where $n$ is even, $h = (b-a)/n$, $x_j = a + jh$. This is equivalent to
% our general weighting scheme where $w_1 = w_n = h/3$, and otherwise $w_j = 4h/3$
% if $j$ is even, and $w_j = 2h/3$ if $j$ is odd. 
%
% All of these rules are fairly easy to program, Simpson's is the most
% common since it turns out to be third-order exact and has an $O(h^4)$
% approximation error. 
%
%
% <include> Int_simp.m </include>
% 
%
% Things are much more accurate: 
Int_simp(@(x) exp(x), 0, 1, 4) - (exp(1) - 1)
Int_simp(@(x) exp(x), 0, 1, 8) - (exp(1) - 1)
Int_simp(@(x) exp(x), 0, 1, 16) - (exp(1) - 1)
Int_simp(@(x) exp(x), 0, 1, 32) - (exp(1) - 1)

%%
% While we can go to higher orders, they typically aren't used. One reason
% is that higher order formulaes will feature negative weights, which is a
% bit strange and can cause numerical problems. 

%% Adaptive Quadrature
% Most implementations of Newton-Cotes methods use adaptave quadrature rules. 
%
% * This is because the
% problem with Newton-Cotes rules is that you provide a number of points to
% use for integral evaluation but that doesn't necessarily give you a good
% idea of the error. 
% * The idea behind these methods is to increase the number of points in
% the calculation and see if the change in the evaluated integral is
% "small"
% * This uses the fact that the error converges to 0 as $N$ rises. 
% * Simpson's rule is very conveinent for this since if you double
% the number of points you can re-use all the previously calculated nodes.
%
% Here's a quick example of Adaptive quadrature: 
%
% <include> Int_Asimp.m </include>
%
% To see how it works: 
[I, N] = Int_Asimp(@(x) exp(x), 0, 10, 1e-8)
I - ( exp(10) - 1)

%%
% * Fancier implementations of adaptive quadrature whill only add points
% where a lot of curvature is detected, rather than sticking to an equally
% spaced grid. 
%
% * MATLAB has a fairly efficient adaptive quadrature method, creatively
% named |integral|.
MatI = integral(@(x) exp(x), 0, 10)
MatI - ( exp(10) - 1)
%
% You can set relative and absolute tolerance for integral, and also
% suggest quadrature nodes, see the documentation.

%% Gaussian Quadrature
%
% Gaussian Quadrature rules are constructed for specific weight functions
% $w(x)$, which often end up being common probability density functions.
% Making them a popular method for computing expectations. 
%
% Rather than "breaking up" the integral, the idea behind quadrature rules
% is to develop a set of nodes and weights which will _exactly_ match all
% polynomials of order-n. The hope is that this will prove to be a good
% approximation for other smooth functions. 
%
% They will be less effective for non-smooth functions, for which
% simulation is an option, or perhaps Newton-Coates methods for functions
% with a small number of kinks or discontinuities (breaking up the function
% will limit the error caused by the kink. 
%
% We can pull it off by solving the following set of moment conditions: 
% 
% $$ \int_I x^k w(x) dx = \sum_{i = 1}^n w_i x_i^k $$
% 
% For $k = 0, \ldots, 2n-1$. So we have $2n$ equations and $2n$ unknowns,
% $(w_i, x_i)_{i=1}^n$. Assuming we can compute the moments of the pdf
% implied by $w(x)$, we just need to solve these non-linear equations. 
%
% For example, suppose we want 3 point quadrature of the weight function of
% the standard normal distribution:
%
% $$w(x) = \frac{1}{\sqrt{2\pi}}e^{-x^2/2}$$
%
% We need to solve this system of equations: 
%
% <include> FindNormPointsWeights.m </include>
%
% Which we can do with a standard nonlinear solver provided we have a
% reasonable starting guess: 
wx = fsolve('FindNormPointsWeights', [.3 .3 .3 -1, 0, 1])

%%
% Now that we have the nodes and weights, we can use them to approximate
% multiple funcitons: 
w = wx(1:3);
x = wx(4:6);

%%
% Suppose I want to approximate a the expected logit choice probability for
% a good of mean value $-1$, but with a heterogeneous taste shock
% distributed Normal(0,1), so 
%
% $$ u_{it} = -1 + x_i + \varepsilon_{it} $$
%
% where $z_i \sim N(0,1)$ and $\varepsilon_{it} \sim T1EV$.
%
% Then my integral is:
% 
% $$ \int \frac{\exp(-1 + x)}{1 + \exp(-1 + x)} \phi(x) dx $$
%
% And using the 3 quadrature points we can approximate it with: 
f = exp(-1+x)./(1 + exp(-1+x))
Ef = w*f'

%% 
% Quadrature weights for most distributions have already been worked out,
% so you can use a package or look them up in a book.  
% The only trick is to make sure you are using the correct change of variables for your problem.
%
% Miranda and Fackler's CE tools has several hand functions for computing
% these: 
addpath('../CEtools/');
[x, w] = qnwnorm(3, 0, 1)
[x, w] = qnwnorm(7, 0, 1)
[x1, w2] = qnwnorm(7, 1, 2)

%%
% Of course, notice the last one is just a linear change of variables in
% the nodes: 
(x*sqrt(2))+1

%%
% They also have functions for lognormal, beta, gamma, and uniform
% distributions. Judd prints several of the quadrature points and weights in his book,
% together with change of variables formulas to help you apply them. 
%
% * Gauss-Chebyshev: Definite integrals on $[-1, 1]$ with weight function $w(x) = (1-x^2)^{1/2}$. Change of variables to unweighted definite integral on $[a,b]$. 
% * Gauss-Legendre: Definite integrals on $[-1, 1]$ with weight function $w(x) = 1$. Change of variables to unweighted definite integral on $[a,b]$. 
% * Gauss-Hermite: Indefinite integrals on $(-\infty, \infty)$ with weight
% funciton $w(x) = e^{-x^2}$. Change of variables to indefinite integrals
% with normal pdf weight function. 
% * Gauss-Laguerre: Indefinite integraols on $[0, \infty)$, with weight
% function $w(x) = e^{-x}$. Directly useful for continuous discounting.
%
% You can find more in the <https://dlmf.nist.gov/3.5#v Digital Library of Mathematical Functions>. 


%% Extension to Multivariate Integrals
%
% Of course, we're likely to have multiple dimensions we want to integrate
% over. Conceptually, the quadrature rules we discussed extend almost
% trivially. In practice, we run into the curse of dimensionality. 
%
% Suppose we want to integrate a function over 2 dimensions, we can just apply single dimensional quadrature to both dimensions:  
%
% $$ \int_{a}^{b} \int_c^d f(x_1,x_2) dx_2 dx_1 \approx \int_{a}^{b} \sum_j w_j f(x_1, x_2^j) dx_1 \approx  \sum_i \sum_j w_i w_j f(x_1^i, x_2^j)$$
%
% So we have the so-called ``tensor product rules''. They would work great
% except that they are very expensive. For example, if you want to use 7
% point quadrature over 5 dimensions you end up having to calculate the
% function $7^6 = 117,649$ times. 
%
% Another option is to try to use <http://sparse-grids.de/ sparse grids>, which are also 
% known as Smolyak points. Heiss and Winschel (2008, _Journal of
% Econometrics_) introduces these methods to economics. I've tried it, but
% wasn't very successful. The issue is that, like higher order polynomial
% approximations some of the weights can be negative. For me, in a 200
% product discrete choice model with 5 dimensions of heterogeneity, 
% I found it was not uncommon for a market
% share approximation to be negative. Maybe there is a creative way around
% this, or maybe it won't be an issue for your application.
%
% The other approach to high dimensionality is to turn to simulation
% methods, which we will discuss next week. 
%
% A second reason to consider simulation is to calculate non-smooth
% integrals. 