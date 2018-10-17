%% Lecture 9 - Monte Carlo Integration
%
% This week, we'll finish numerical integration with Monte Carlo methods.
% And also review a quick application where we use numerical integration to
% estimate a bivariate probit. 
%
%% Monte Carlo Integration: Basic Idea
% Use these when: 
%
% * Function is high dimensional. 
% * Function is non-smooth or discontinuous (indicator funcitons). 
%
% Esentially, we are just leveraging the law of large numbers. If $x_s$ are
% iid realizations of the random variable $X$ and $f$ is "nice",
%
% $$ \lim_{S \rightarrow \infty} \frac{1}{S} \sum_{s = 1}^S f(x_s) =
% Ef(X) $$
% 
% So monte carlo integration just amounts to drawing a random sample from
% the distribution of interest and taking an average. 
%
% There are just two things do deal with:
%
% * How do you create random draws.
% * How do you draw from the correct distribution.
%
%% Generating "Random" Numbers
%
% A computer is a deterministic machine, so it cannot really produce random
% numbers. Most commonly it uses psueo-random numbers based on
% deterministic sequences which ``look'' random:
%
% * Zero serial correllation (at all lags)
% * Correct frequency of "runs" (e.g., three draws greater than .5). 
%
% Probably the simplist of these sequences is the <https://en.wikipedia.org/wiki/Linear_congruential_generator Linear Congruential
% Generator> (LCG). Which is a sequence of the form: 
%
% $$ X_{n+1} = (aX_n + c) \mbox{ mod } m $$
%
% * $a$, where $0 < a < m$, is the multiplier.
% * $c$, where $0 < c < m$ is the increment. 
% * $m$, is the modulous.
% * $X_0$, is the initial seed.
%
% These simple RNGs are no longer commonly used. In general, you should use
% the random number generators available in MATLAB or C packages. Since you
% can set different seeds, there are plenty to choose from. 
%
% Just for fun, here is a toy LCG-RNG:
a = 65539;
m = 2^31;
c = 1;

seed = 2012;
Seq = zeros(5000, 1);
Seq(1) = seed;
for s = 1:5000-1
    Seq(s+1) = mod(a*Seq(s) + c, m);
end

%Now convert to [0,1]:
Seq = Seq./(m-1);

mean(Seq)

%%
%Let's split it into two dimensions and plot it:
X1 = Seq(1:2:end);
X2 = Seq(2:2:end);
scatter(X1, X2);

%% 
% Looks pretty random to me. And in fact this is very close to the <https://en.wikipedia.org/wiki/RANDU RANDU> LCG commonly
% used in the 1960 and 70s, which turns out to be one of the *"worst RNGs
% ever deisgned".* I've helped it out slightly by setting c=1 instead of
% c=0, but it's still not something I would recommend even though it "looks"
% random from this plot. 
% 
% The main lesson is not to code these things yourself, since the accepted
% versions have been run through many tests for serial correlation that you
% won't have time or interest in conducting. 
%
% Instead, just use your local package. Since we are using MATLAB, the
% native random number generator works just fine. For most applications, you can just us rng to
% set the state and your preferred generator (just use the default, Mersenne Twister).
%
% Set the seed: 
seed = 8673310; 
rng(seed); 
Y = rand(2,2)
Z = rand(2, 2)

%% 
%
% Setting and saving the seed lets me "reset" the generator to its initial state so I get the
% same "random" sequence again. This is important so that my results are
% reproduceable. 
rng(seed);
rand(2,2)

%% Non-Uniform Distributions
% So what if you want to draw from something other than a Uniform
% distribution? 
%
% In general, you can handle this using cumulative distributions functions. The CDF of $X$,
%
% $$ F(x) = pr(X \leq x), $$
%
% maps a random variable to a Uniform(0,1). Therefore if we have a uniform
% draw, the inverse CDF maps it to a draw distributed according to F. So if
% $u \sim Unif(0,1)$ then
%
% $$x = F^{-1}(u) $$ 
% 
% is distributed according to $F$. If $F$ needs to be approximated, we can
% do that, although things get harder if $F$ does not have a closed form and is multi-dimensional.
% For that we'll talk about MCMC methods later in the course. 

% So if we want a normal random variable, we can just draw a uniform run
% it thorugh the inverse normal cdf:
normY = norminv(Y)
%%
% But we can also just draw standard normal random variables directly: 
normZ = randn(2,2)

%% 
% Frequently, we will want to draw from a jointly normal distribution, so
% it is useful to know how to construct these from independent random
% normals. This is just a matrix multiplication. 
%
% Suppose a variable is distributed $N(\mu, \Sigma)$, where 
%
% $$ \mu = \left( \begin{array}{c} \mu_1 \\ \mu_2 \end{array} \right) $$
%
% $$ \Sigma = \left( \begin{array}{cc} \sigma_1^2 & \rho \sigma_1 \sigma_2 \\
%                                   \rho \sigma_1 \sigma_2 & \sigma_2^2 \end{array} \right) $$
%
% The Cholesky factorization of $\Sigma$ decomposes it into a triangular
% matrix and its transpose, $\Sigma = LL'$,
Mu = [2; 3];
Sigma = [ 1 .9; .9 2];
U = chol(Sigma) %Annoyingly, MATLAB returns the upper matrix by default, just be aware. 
U'*U

%%
% Now we can draw 500 bivariate normal by first drawing standard normals and
% then using $L$: 
Z = randn(2, 2500);
biNorm = repmat(Mu, 1, 2500) + U' * Z;
scatter(biNorm(1,:), biNorm(2,:));

%% 
% It is also possible to directly parameterize the Cholesky
% decompositition, which is convenient if you need to derive the derivatives of your estimator.  
% Just be aware when you do so that those parameters are
% not interpretable as variances and covariances. 
%
% $$ \left( \begin{array}{cc} a &  0 \\
%                                   b & c \end{array} \right)
%  \left( \begin{array}{cc} a &  b \\
%                                   0 & c \end{array} \right) = 
%  \left( \begin{array}{cc} a^2 &  ab \\
%                                   ab & b^2 + c^2 \end{array} \right) = 
% \left( \begin{array}{cc} \sigma_1^2 & \rho \sigma_1 \sigma_2 \\
%                                   \rho \sigma_1 \sigma_2 & \sigma_2^2 \end{array} \right) $$
%
% Since we have three equations and three unknowns, we can solve for
% $\sigma_1, \sigma_2, \rho$ from $a, b, c$, of course, we'll need to use
% the delta method to do inference. Some people just report the Cholesky
% matrix itself. 
%
%% Quasi-Monte Carlo
%
% While pseudo-Monte Carlo relies on the Law of large numbers, with
% inference based on the central limit theorem. It can be slow, since a CLT
% exhibits $\sqrt{n}$ convergence, we need to quadruple the sample to halve
% the error. 
%
% Quasi-Monte Carlo sequences are deterministic, but concentrate on
% efficiently, they effectively dispense with the effort to mimic random
% numbers and instead concentrate on what we want from integration, that sequences
% are \emph{equidistributed}. A sequence $\{x_j\} is equidistributed over [a, b] iff::
%
% $$ \lim_{n \rightarrow \infty} \frac{b - a}{n} \sum_{j=1}^n f(x_j) = \int_a^b f(x) dx $$
% 
% Which doesn't say anything about the serial correlation between draws.
% Effectively, Monte Carlo methods rely on the Law of large numbers, which
% says that random sequences will be equi-distributed. Quasi-monte carlo
% gives up on randomness and uses number theory to find equidistributed
% deterministic sequences. (In actually has nothing to do with randomness,
% or Monte Carlo, at all). 
%
% There are only a handful of known deterministic sequences that are proven to be equidistibuted.
% Fortunately, that's all we really need. Monte Carlo evidence has shown
% that for the same number of points quasi-MC sequences can be up to 10x
% more accurate than pseudo-MC sequences. This makes them the most
% efficient thing going when Gaussian quadrature isn't an option. 
%
% The formulas for Weyl, Harber, Niederreiter, and Baker sequences are
% provided in Judd. <https://en.wikipedia.org/wiki/Halton_sequence Halton sequences>  is another variant which I've used.
% Sobol points are also popular. 
%
% All of these points are relatively simple to program. But there are many
% packages available to produce them. 
%
% The CEToolbox contains an algorithm to produce `N'eiderrieter, 'W'eyl and 'H'aber points
% Do to draw 500 Neiderrieter points on the unit square:
addpath('../CETools');
[n, w] = qnwequi(500, [0 0], [1, 1], 'N');
[n(1:5, :) w(1:5)]
scatter(n(:,1), n(:,2));

%%
% I'm also providing some MATLAB code for Halton Sequences, this code
% doesn't provides weights (since they are just 1/n) and only returns
% values for Uniform [0,1], since you can do the transformation to [a, b]
% yourself. 
h = haltonseq(500,2);
scatter(h(:,1), h(:,2), 'k');

%%
% Compare to psuedo-random draws:
y = rand(500,2);
scatter(y(:,1), y(:,2), 'r');

%%
% Of course, we usually want to integrate over something other than uniform
% distributions. No problem, we can use the cdf-inverse to transform the
% points just as though they were pseudo-random points. I've provided
% |HaltonNormShuffle| which does this for Halton sequences. It also
% "shuffles" the points to reduce colinearity across dimensions in Halton
% sequences. 

hNorm = haltonNormShuffle(500, 2, 33456); %Last argument is a seed for shuffling dimensions. 
scatter(hNorm(1,:), hNorm(2,:), 'k');
%%
% We can convert these from $N(0,1)$ to $N(\mu, \Sigma)$ just as we did for
% pseudo-MC draws. 

hbiNorm = repmat(Mu,1,500) + U' * hNorm;
scatter(hbiNorm(1,:), hbiNorm(2,:));


%% Using Numerical Integration in Estimation
%
% Let's look at an application from Roberts and Tybout (1997, AER). They
% have data from 650 manufacturing plants in Colombia and want to estimate
% a model of plant exporting that accounts for persistence in exporting
% behavior and unobserved plant heterogeneity. 
%
% For the years 1984-1989 the long run profit from exporting is modeled as a
% function of plant characteristics and lagged values of the export dummy.
% The idea is that long run profits from exporting can be modeled as,
% 
% $$ a_{it}=\mu _{0}+\gamma _{1}Y_{it-1}+\gamma_{2}{Y}_{it-2}+\gamma_{3}{Y}_{it-3}+\sum_{t=85}^{89}\mu_{t}D_{t}+\beta Z_{it}+\alpha _{i} $$
%
% Where $Y_{it}$ is an export summy, $D_t$ are year dummies $Z_{it}$ are
% observable plant charachteristics, and $\alpha_i$ is a plant effect that is constant over time.
% 
% If we add a standard normal shock to $a_it$ to determine the export
% decision we get that firm $i$ will export in year $t$ with probability: 
%
% $$ \mbox{Pr}(Y_{it} = 1) = \mbox{Pr}(a_{it} + \varepsilon_{it}  \geq 0) =
% \Phi(a_{it}) $$
%
% However, this model suffers from an initial conditions problem. In the
% first year (1984) the lagged values will be correllated with the plant
% effect $\alpha_i$, but this is not modeled. We'll deal with this with a
% correction proposed by Heckman (1981). For years the initial condition
% years we'll propose a seperate model for long run profits that does not
% include lagged years (which we don't observe). 
%
% $$ b_{it}=\lambda _{0}+\sum_{t=82}^{83}\lambda _{t}D_{t}+\lambda
%  Z_{it}^{P}+\alpha _{i}^{P} $$
%
% As above, we'll assume a normal chock to long-run profits determines the
% exporting decision, $\mbox{Pr}(Y_{it} = 1) = \Phi(b_{it})$.
% There are no lagged dependent variables here so no initial conditions problem. 
% Moreover, we will allow
% $\alpha_i$ and $\alpha_i^P$ to be correlated, which will control for the correlation in the
% equation for $a_{it}$. In fact, we'll assume they
% are jointly normal with variance matrix: 
%
% $$ \Omega = \left( \begin{array}{cc} \omega_{11} &  \omega_{12} \\
%                                   \omega_{12} & \omega_{22} \end{array}
%                                   \right) $$
%
% To build the likelihood, start with the likelihood of a single plant $i$,
% conditional on $(\alpha_i^P, \alpha_i)$, 
%
% $$ L(\mu, \gamma, \beta, \lambda| \alpha_i^P, \alpha_i ) = \left(
%  \Pi _{t=81}^{83}\left[ \Phi (b_{it})^{Y_{it}}(1-\Phi (b_{it}))^{(1-Y_{it})}
%  \right] \right)  \cdot \left( \Pi _{t=84}^{89}\left[ \Phi (a_{it})^{Y_{it}}(1-\Phi
%  (a_{it}))^{(1-Y_{it})}\right] \right). $$
%
% Of course, we don't observe $(\alpha_i^P, \alpha_i)$, so they will have
% to be integrated out. We assume shocks are independent across plants so the full data likelihood is, 
%
% $$ L(\mu ,\gamma ,\beta ,\lambda ,\Omega )=\Pi _{i=1}^{650}\int_{-\infty}^{\infty }\int_{-\infty }^{\infty }L_{i}(\mu ,\gamma ,\beta ,\lambda|\alpha _{i},\alpha _{i}^{P})
%   \phi (\alpha ,\alpha ^{P}|\Omega )d\alpha d\alpha ^{P} $$
%
% So we need to use numerical integration to optimize this. It seems like a good candidate for Gaussian Quadrature, so let's use that. Here is the likelihood function: 
%
% <include> llk.m </include>
%
% Some things I might make an assignment out of: 
%
% * Notice that the code parameterizes the Cholesky root, then creates the
% varience-covariance to pass to qnwnorm. Mildly inefficient. 
% * Might be a good idea to avoid calling qnwnorm every funciton eval.
% Because this is optimized via Nelder-Mead and requires a lot of
% evaluations, every bit matters. 
% * If we do the two above, we could easily alter this code to use other integration schemes. How? 
% * I've included two files to optimize this, |ex6mod.m| just calls
% Nelder-Mead once, |ex6.m| loops. Why?
% * If you play with this, don't forget to add the |CETools| library to
% your path. 

