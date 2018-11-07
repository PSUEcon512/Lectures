%% Lectures 10 & 11: Differentiated Products Demand with Random Coefficients (BLP)
%
% The differentiated products demand system is a workhorse model for
% product competition in IO. Developed in Berry, Levinsohn, and Pakes (1995, Ema), 
% We'll use the "practitioners guide" paper by Nevo (2000, JEMS) as our
% guide. I'll go over code which was developed by Chirs Conlon at NYU,
% which is updated from Nevo's MATLAB code. 
%
% While this model is useful to know in and of itself it also features all
% the elements we have covered so far this semester:
% 
% * Solving linear and non-linear equations. 
% * Optimization.
% * Numerical integration.
%

%% Data & Setup
%
% The aim of this model is to understand market-level substitution patterns
% through the lends of a discrete choice model with heterogeneous agents.
% We observe: 
%
% * $t = 1,\ldots,T$ *markets*, each containing a mass of consumers. For
% asymptotics, we will assume $T$ grows to infinity. These may be repeated
% observation of the same market over time (though we assume the model is
% static) or data from distinct markets where the same products compete
% (e.g., different MSAs in the United States). 
% * A set of $J$ *products* which are offered in the markets. Each product
% has an average price in the market $p_{jt}$ and charachteristics which
% may vary by market $x_{jt}$. Not all products need to be offered in every
% market. (In fact, it may be helpful if there is variation in the choice
% set, assuming it is exogenous. 
%
% We assume that consumers buy at most one product, and they do so to
% maximize their indirect utility conditional on purchase, where indirect
% utility of consumer $i$ in market $t$ of product $j$ is: 
%
% $$ u_{ijt} = \alpha_i(y_i - p_{jt}) + x_{jt} \beta_i + \xi_{jt} +
% \varepsilon_{ijt} $$
% 
% * $y_i$ is consumer i's income, so $(y_i - p_{jt})$ represents consumers
% i's expenditure on other goods (outside this market). The coefficient
% $\alpha_i$ is consumer i's price sensitivity (so allows for
% non-homothetic preferences). 
% * $\beta_i$ represents consumer i's taste for observable characteristic
% $x_{jt}$. 
% * $\xi_{jt}$ is an unobserved quality of product $j$ in market $t$ that
% is common to all consumers in $j$. It is important as it is potentially
% correlated with price. Importantly, it will be assumed to be
% uncorrellated with $x_{jt}$.
% * $\varepsilon_{ijt}$ is a consumer-specific taste shock. Since it is
% specific to a consumer it is reasonable to assume that it is uncorellated
% with product charachteristics. 
%
% We aren't going to observe individual level data, only market level data,
% however we may assume that we learn the _distribution_ of consumer level
% charachteristics across markets. E.g., we can learn the income and age
% distribution by MSA. Therefore, we will assume that consumer tastes can
% be modeled according to distributional assumptions: 
%
% $$ \left( \begin{array}{c} \alpha_i \\ \beta_i \end{array} \right) = 
%     \left( \begin{array}{c} \alpha \\ \beta \end{array} \right) + \Pi D_i
%     + \Sigma \nu_i $$
%
% * $D_i \sim P_{Dt}(D)$ represents a draw from the distribution of
% demographics (income, race, family size) it is a $Z x 1$ vector. 
% * $\Pi$ is a $KxZ$ matrix that determines how demographic characteristics
% affect tastes on average (e.g., families with kids like sugary cereals
% more). 
% * $\nu_i \sim P_\nu(\nu)$ represents a draw from the distribution of
% unobserved tastes. It is a $K x 1$ vector. The distribution $P_\nu$ is
% assumed to be known (typicaly standard normal or log-normal). 
% * $\Sigma$ is a KxK vector that represents the Cholesky decomposition of
% the variance-covariance matrix of unobserved tastes. This has important
% implications for the estimation of substitution patterns. If $\sigma_kk$
% is very large, it means products with a lot of charachteristic $k$ will
% be unlikely to subsitute to products with little $k$. Very often,
% $\Sigma$ is a diagonal matrix to reduce the number of parameters to
% estimate. 
%
% In addition to all the products that are purchased, we consider an
% outside good or ``no purchase'' option. This will simply be considered
% good 0: 
%
% $$ u_{i0t} = \alpha_i y_i + \xi_{0t} + \pi_0D_i + \sigma_0 v_{i0} +
% \varepsilon_{i0t} $$
%
% * The parameters $\pi_0$ and $\sigma_0$ allows for the outside good to have its own demographic and random
% effects however they are not identified seperately from an constant term
% on the inside goods, so we usually just set them to 0. (Not a huge deal
% since these are essentially the same thing). 
% * Similarly, a $\xi_{0t}$ is not identified unless we normalize one of
% the inside goods across markets. 
%
% Given this structure, we can re-write indirect utilities in terms of
% linear and nonlinear parameters. Where linear parameters impact the mean
% utility of a product in market $j$ and non-linear parameters affect
% deviations from that mean. 
%
% $$ u_{ijt} = \alpha_i y_i + \delta_jt(x_{jt}, p_{jt}, \xi_{jt}; \theta_1)
%    + \mu_{ijt}(x_{jt}, p_{jt}, \nu_{i}, D_i; \theta_2) + \varepsilon_{ijt} $$  
%
% * $\delta_{jt} = x_{jt} \beta - \alpha p_{jt} + \xi_{jt}$ is the utility of the consumer with "mean 0" demographics and taste shocks.
%   $\Theta_1 = (\alpha, \beta)$, the mean tastes. 
% * $\mu_{ijt} = [-p_{jt}, x_{jt}] ( \Pi D_i + \Sigma v_i )$ is this
% consumers i's deviation in utility for product j from mean, it is a function of $\theta_2 = (\Pi, \Sigma)$.  
%
% Why are we doing this? Because we will exploit the result that, for a
% given distribution of heterogeneity, there will be a one to one mappint
% from observed shares $s_{\cdot t}$ to mean utilities in the market
% $\delta_{\cdot t}$. That will allow us to create moments that are
% additively seperable in the structural error $\xi_{jt}$. 
%
%
% We observe market shares of each product in each market. To determine
% these, we need to integrate over all consumer-level charachteristics.
% Specifially define A_{jt} as the set of consumers who buy product $j$:
%
% $$ A_{jt}(x_{\cdot t}, p_{\cdot t}, \delta_{\cdot t}; \theta_2) = \left
% \{ (D_i, \nu_i, \varepsilon_{i, \cdot, t}) : \forall \ell \in J, u_{ijt}
% \geq u_{i \ell t} \right\} $$
%
% * Notice we right this just as a function of delta and deviations from
% mean utility. That is one reason for this notation. 
% * Now to calculate the predicted market shares, we just need to integrate over the distribution of tastes: 
%
% $$ s_{jt} = \int_{A_{jt}} dP_\varepsilon(\varepsilon) dP_\nu(\nu) dP_D(D) $$
%

%% Distributional Assumptions & Substitution Matrix
%
% * We will assume that $\varepsilon$ is a Type-I Extreme value shock, so
% it can be integrated out analytically. 
%
% * Demographics are typically either simulated from a non-parametric
% distribution or from a parametric fit from data. 
%
% * Taste shocks $\nu_i$ will be numerically integrated from an assumed
% standardized distribution. 
%
% The main thing this work buys us is flexible (and hopefully realistic)
% substitution patterns. What do we mean by this? Suppose we dropped
% demographics and taste heterogeneity from the model. So the only thing
% allowing consumers to make different choices was $\varepsilon_{ijt}$.
% Then we'd have the classic multinomial logit model: 
%
% $$ s_{jt} = \frac{ \exp ( x_{jt} \beta - \alpha p_{jt} + \xi_{jt} ) }{ 1 +
% \sum_{\ell \in J} \exp ( x_{\ell t} \beta - \alpha p_{\ell t} + \xi_{\ell
% t} ) } $$
%
% This involved no numerical integration, but also has the independence of
% irrelevant alternatives property, which implies substitution pattersn are essentially a function of shares. 
% In fact the formula for elasticities are: 
%
% $$ \eta_{jkt} = \frac{\partial s_{jt} p_{kt}}{\partial p_{kt} s_{jt}} =  \left\{ \begin{array}{lr} -\alpha p_{jt}(1-s_{jt}) & \mbox{if } j = k \\ 
%                                  \alpha p_{kt} s_{kt} & \mbox{if } j \neq k  \end{array} \right. $$
%
% So all the work of determining price elasticities is done by a single
% parameter, $\alpha$, for both own and cross-price substitution. However,
% in the full model, we have a logit-type substitution pattern at the
% consumer level, but must integrarte over consumers to get the market
% level elasticities the firms care about: 
%
% $$ \eta_{jkt} = \frac{\partial s_{jt} p_{kt}}{\partial p_{kt} s_{jt}} =  \left\{ \begin{array}{lr} - \frac{p_{jt}}{s_{jt}} \int \alpha_i s_{ijt}(1-s_{ijt}) dP_D(D) dP_\nu(\nu) & \mbox{if } j = k \\ 
%                   \frac{p_{kt}}{s_{jt}} \int \alpha_i s_{ijt} s_{ikt} dP_D(D) dP_\nu(\nu) & \mbox{if } j \neq k  \end{array} \right. $$
%
%% Overview of Estimating the Model 
%
% The first thing one might try is a simple non-linear least squares to match observed shares: 
%
% $$\min_\theta || s(x, p, \delta(x, p, \xi,; \theta_1); \theta_2) - S|| $$
%
% * First big problem is $\xi$, we haven't put a distributional assumption
% on that though we could. 
% * But then we'd run into the issue that $p$ is in principle correlated
% with $\xi$. How? We'd need a supply side and more assumptions. 
% * Also ALL of the parameters enter this optimization non-linearly, as
% opposed to just $\theta_2$ in the method we will use. 
%
% Instead, we use the following approach. We assume we have instruments
% which are uncorrelated with the structural error such that 
%
% $$ E[Z\xi(\theta^*)] = 0 $$
%
% Where $\xi(\theta^*)$ is the implied structural error computed for the
% true value of $\theta$. (I.e., $\xi$ is the true structural error.)  The
% important thing is this is a linear moment condition, thought
% $xi(\theta^*)$ must be computed nonlinearly. 
%
% This leads to the following method:
%
% # Fix $\theta_2$, there is a 1-to-1 mapping from shares to mean
% utilities, so solve the non-linear equation $S - s(x, p, \delta;
% \theta_2)$ for $\delta$.
% # Now we have the equation, $\delta_{jt}(\theta_2) = x_{jt}\beta - \alpha p_jt +
% \xi_{jt}$, a linear IV regression can be used to determine
% $\xi_{jt}(\theta_2)$. 
% # Search over $\theta_2$ to minimize the GMM objective: 
%
% $$ \min_{\theta_2} \xi(\theta_2)'Z\Omega Z' \xi(\theta_2) $$
%
%
%
%% Computing Market Shares:
%
% We first need a way to compute $s(x, p, \delta; \theta_2)$. This is
% essentially an applciation of numerical integration. 
%
% Let $S$ be the number of nodes we draw from the joint distribution of
% tastes and heterogeneity. That is we draw $\left\{ (\nu_i, D_i)
% \right\}_{i=1}^S$. If we are using quadrature, we have weights $w_i$. 
%
% $$ s(x, p, \delta, P_{ns}; \theta_2) = \sum_{i = 1}^S w_i
% \frac{\exp \left\{ \delta_{jt} + [-p_{jt}, x_{jt}](\Pi D_i + \Sigma v_i) \right\} }
% { 1 + \sum_{\ell \in J_t} \exp \left\{ \delta_{\ell t} + [-p_{\ell t}, x_{\ell t}](\Pi D_i + \Sigma v_i) \right\}} $$
%
% This is a straightforward calculation, but we do have to think about
% overflow/underflow issues since we are exponentiating. Overflow is the
% more likely problem, and can be dealt with by just normalizing by the
% maximum indirect utility. 
%
% $$ \frac{ \exp(a + b)}{1 + \exp(a + b)} = \frac{ \exp(a)}{\exp(-b) +
% \exp(a)} $$
%
% So we can usually keep all the utilities in the computable range. Here is
% a look at conlon's share calculation, utlity normalization is actually
% commented out. 
%
% <include> blp-conlon/rc_share_safe.m </include>
%

%% Inverting shares to determine mean utilities: 
%
% Next, we need to be able to solve the inversion of shares to mean
% utilities. So we have find $\delta$ to solve the nonlinear equation
%
% $$ S - s(x, p, \delta; \theta_2) = 0 $$
%
% This must be done once for each candidate parameter vector during
% estimation. On the bright side, it can be parallelized to be computed
% market by market. 
%
% BLP note that this problem can be solved using using a contraction:
% 
% $$ \delta_{\cdot t}^{h+1} = \delta_{\cdot t}^{h} + \log S_{\cdot t} - \log
% s_{\cdot t}(x_{\cdot t}, p_{\cdot t}, \delta_{\cdot t}^h; \theta_2) $$
%
% If we have a good start point for $\delta$, it can also be solved using
% Newton's method. 
%
% People have been staring at this problem for some time. Recently
% Reynaerts, Varadhan and Nash (2012), proposed an acceleration based on a
% squared polynomial extrapolation of the draws (SQUAREM)
% It remains globally convergent (an issue with Newton's method 
% while being 5 times faster than the "traditional" contraction. It is
% Conlon's preferred implementation. 
%
% This code shows parallelization by market, and the various options for
% solving the fixed point: 
%
% <include> blp-conlon/solveAllShares.m </include>
%
% And this function shows the two contraction methods:
% |blp-conlon/fp_squarem.m|. It is written as a generic contraction
% accelerator. Hence we pass the BLP contraction as an argument.
%
% Notice that this code also computes the Jacobean of the mean utility
% vector with respect to parameters, we'll see later where this is useful.

%% IV regression of mean utilities to recover unobserved quality:
%
% Now for a given $\theta_2$ we know how to determine $\delta(\theta_2)$,
% the beauty of this is that now we can recover $\xi(\theta_2)$ as the
% residual of a linear IV regression: 
%
% $$ \delta_{jt}(\theta_2) = x_{jt} \beta - \alpha p_{jt} + \xi_{jt}(\theta_2) $$
%
% Where the instrument set $z_{jt} = (x_{jt}, \tilde{z}_{jt})$ where
% $\tilde{z}_{jt}$ is a vector of valid instruments.  Note that we need
% enough instruments to identify ALL parameters in the model:
%
% * The $x_{jt}$ instruments provide moments to identify the $\beta$
% coefficients. 
% * One valid instrument for price will identify $\alpha$. 
% * If this IV regression is just identified, our GMM objective function
% will be 0, what is left to identify $\theta_2 = (\Pi, \Sigma)$? Need
% additional instruments to identify these parameters. 
%
% Since this is a computational course, we'll just assume such instruments
% exist. Finding them is a big part of doing reasonable empirical work. 
%
% Conlon's IV regression is quite straightforward, just need to call |[beta,resid]=ivregression(delta,X,Z,W);|: 
%
% <include> blp-conlon/ivregression.m </include>
%
% You could save yourself solving a linear equation every iteration by
% computing the projection matrix and turning this into a matrix
% multiplication, but this isn't the time sink of the algorithm (solving
% for delta is) and this is a bit more intuitive. 
%

%% Compute GMM Objective:
%
% Now we've got everything we need to compute the GMM Objective function,
% which is just a matrix multiplicaiton. 


%%
%   function [fval,g]=evalSingle(theta)
%     % this extracts the parameters for your specification
%     p = get_params(theta,draws);
%     [delta,Jac]=solveAllShares(dtable,draws,p,method);
%     dtable.delta=delta;
%     [beta,resid]=ivregression(delta,X,Z,W);
%     fval=(resid'*Z)*W*(resid'*Z)';
%     g=-2*(Jac'*Z)*W*(resid'*Z)';
%   end


%%
% Where are |X|, |Z|, |W|, |draws|, and |dtable| coming from? 
%
% This code is using nested functions, so when we look at
% |solveRCBLPpar(dtable,draws,theta0,extract_fun)| you will see that |evalSingle|
% is defined inside, and so  inherits its scope. It is
% a nice trick to keep things clean. 

%% Optimization and Two-Step GMM
%
% Once we have the objective function (and its Jacobian), we just need to
% use an optimizer to solve it. 
%
% * The code is set up to use either |fmincon|
% or |knitromatlab|. 
% * I have removed |knitromatlab| since the podium machine
% doesn't have license access.
% * Notice: As we discussed in the theory portion, optimization is only
% over nonlinear parameters, so we use those to recover the nonlinear
% parameters post-Optimization. 


%%    
%   function [res]=get_results(tableA,x0)
%    % function handle f is mapped to evalsingle below for a (X,Z,W) 
%    if 0 %KNTRO license currently down, %exist('knitromatlab'),
%        [that]=knitromatlab(f,x0,[],[],[],[],lb,ub,[],[],ops);
%     else,
%        [that]=fmincon(f,x0,[],[],[],[],lb,ub,[],ops);
%     end
%     %
%     % After optimization recover the linear parameters and objective 
%     thetahat =get_params(that,draws);
%     delta=solveAllShares(tableA,draws,thetahat,method);
%     dtable.delta=delta;
%     [beta,resid]=ivregression(delta,X,Z,W);
%     fval=(resid'*Z)*W*(resid'*Z)';
%     %
%     % Put the results into structure
%     [~,SEest]=getCovariance(that,dtable,draws);
%     res.fval = fval; res.beta=beta; res.theta = that; 
%     res.delta=delta; res.resid = resid; res.SE=full(SEest);
%   end

%%
% Finally, since we have an overidentified GMM problem, we need to use the
% optimal weight matrix. Once we have a pilot estimate we can construct the
% optimal weight matrix
%
% $$ \hat{W}_{opt} = (Z'\hat{\xi}\hat{\xi}'Z)^{-1} $$
%
% and optimize a second time. Any weight matrix is consistnent for the
% first stage, in practice we use $W_0 = (Z'Z)^{-1}$.
%
% The two-step procedure is the heart of |solveRCBLPpar| (par for parallel,
% though I have disabled that portion of the code). It is quite
% anti-climactic:

%%
%   tic
%   % first step
%   [results1]=get_results(dtable,theta0);
%   print_results(results1);
%   % update weight matrix and produce second-step
%   [W,~]=getCovariance(results1.theta,dtable,draws);
%   [results2]=get_results(dtable,results1.theta);
%   print_results(results2);
%   toc

%%
% And that's it! Now we are ready to take a stroll through the code. 