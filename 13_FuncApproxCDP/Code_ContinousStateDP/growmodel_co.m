%Let's do value function iteration with collocation on a simple capital
%growth model. 

% Section 1: Functional Forms
clear;

% the model is the simple growth model  
% preferences: u(c)=c to the power of (1-sigma) over (1-sigma)
% beta is the discount factor
% delta is the rate of capital depreciation
% the production function is: k to the power of alpha

% Section 2
% Basic Parameters

     m.sigma=1;
     m.beta=.9;
     m.alpha=.75;
     m.delta=0.3;
     
%% Section 3: Define Spaces and Initialize

% determine the steady state capital stock 
% so here Kstar is something we solved for analytically

Kstar=(((1/(m.alpha*m.beta))-((1-m.delta)/m.alpha)))^(1/(m.alpha-1));

% put in bounds here for the capital stock space
% We'll be using Chebyshev Interpolation between these bounds.
cheb.lo=Kstar*.8;
cheb.hi=Kstar*1.2;
cheb.order = 10;
[Knodes, Phi] = setupCheb(cheb.order, cheb.lo, cheb.hi);

%
vcur = zeros(size(Knodes));
vnew = vcur;
vmax = vcur;

%Initial value function is all zeros
c = zeros(cheb.order, 1);

%Test of veval
%veval(Knodes(1), Knodes(2), c, cheb, m);

%% Section 4: Value Function Iteration with Collocation
maxit = 200;
options = optimset('Display', 'off', 'Algorithm', 'sqp');
tic;
for iter = 1:maxit

    %For each node, optimize
    for kk = 1:length(Knodes)
       [kp(kk) vmax(kk)] = fmincon(@(k) veval(k, Knodes(kk), c, cheb, m),Kstar, [],[], [], [], 0, Knodes(kk).^m.alpha + (1-m.delta)*Knodes(kk), [], options);
    end
    
    %Now use new vmax to form next approximation of the value function:
    vnew = -vmax;
    cnew = Phi\vnew;

    %Now Check and see if we've converge, if not, continue loop.
    tol = norm(vnew - vcur);
    if tol<1e-5
        disp(sprintf('Convergence Achieved, %d Iterations', iter));
        toc;
        break
    else
        if mod(iter, 10)==0
           disp(sprintf('Iter: %d, tol = %f', iter, tol));
           disp(c)
        end
        c = cnew;
        vcur = vnew;
    end
end

%% Section 5: Plots

close all;
%Now, Plot the value function:
figure(1);
scatter(Knodes, vnew)
hold on;
K = [cheb.lo:.05:cheb.hi]';
V = chebEval(K,cnew, cheb.order, cheb.lo, cheb.hi);
plot(K,V);
title('Value Function');
xlabel('current capital');

%Next, Plot the Policy Function:
kp_c = Phi\kp';
figure(2);
scatter(Knodes, kp);
hold on;
Kp = chebEval(K,kp_c, cheb.order, cheb.lo, cheb.hi);
plot(K,Kp)
plot(K,K,'r');
xlabel('current capital');
ylabel('next period capital');
title('Policy Function');
