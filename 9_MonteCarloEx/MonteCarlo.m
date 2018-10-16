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
% realizations of the random variable $X$ and $f$ is "nice", 
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
%The seed lets me "reset" the generator to its initial state so I get the
%same random sequence again. This is important so that my results are
%reproduceable. 
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
% Suppose a variable is distributed $N(\mu, \Sigma)$.  

%% Quasi-Monte Carlo

%% Using Numerical Integration in Estimation