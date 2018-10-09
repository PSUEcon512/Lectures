%% Lecture 8 - Numerical Integration
%
% For the next two weeks, we'll consider how to evaluaate the definite integral of a
% real valued function over some interval $I \in R^d$:
%
% $$ \int_{x \in I} f(x) w(x) dx $$
%
% The weight function $w(x)$ will most commonly be a probability density function,
% in which case the integral denotes the expectation of $f(x)$. 
%
% There are plenty of places where such a problem will pop up:
%
% * Evaluating expectations over states of the world.
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
% We can go further by approximating the function with piecewise quadratic
% funcitons instead of piecewise linear, which leads to *Simpson's Rule*. To fit a quadratic on the interval
% we need three points, so we bring back the midpont: 
%
% $$ \int_a^b f(x) dx = \frac{(b - a)}{6} [f(a) + 4f\left(\frac{a+b}{2} \right) + f(b)] - 
% \frac{(b-a)^5}{2880}f''''(\xi) $$
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
% While we can go to higher orders, they typically aren't used. Instead,
% it is more common to use *adaptive quadrature*. This typially raises the
% number of points in the calculation until the integral approximation
% converges. Simpson's rule is very conveinent for this since if you double
% the number of points you can re-use all the previously calculated nodes.
%

%% Change of Variables for Infinite Domains

%% Gaussian Quadrature



