% This is the main file for the estimation of the wind-border paper by
% Cosar, Grieco, and Tintelnot
diary
clear;
clc;

%Either run these scripts or read in their output...
%setup;
load data_structures


%Define initial guesses

rho_ger = ones(m.num_pro_ger * (m.num_firms_ger + 1) ,1) ./ (m.num_firms_ger + 1);    
rho_dnk = ones(m.num_pro_dnk * (m.num_firms_dnk + 1) ,1) ./ (m.num_firms_dnk + 1);    
fe      = ones(m.num_firms_ger,1);         % firm fixed effect (note: all the firms are active in Germany)
betaLogD = .1;
betaFor  = .7;
%betaState = 1.4;


x_0 = [rho_ger;rho_dnk;fe;betaLogD;betaFor];

%% Sparsity pattern

%These are the sparsity templates for each project the first row, is the 
%fringe share and the final column is theadding up constraint
J_m_ger = [[ones(1,m.num_firms_ger);eye(m.num_firms_ger)],ones(m.num_firms_ger+1,1)];
J_m_dnk = [[ones(1,m.num_firms_dnk);eye(m.num_firms_dnk)],ones(m.num_firms_dnk+1,1)];

%This is the Jacobian of the Constraints by the market share parameters
%(rhos)
JJ = [kron(eye(m.num_pro_ger),J_m_ger), zeros((m.num_firms_ger+1)*m.num_pro_ger, (m.num_firms_dnk+1)*m.num_pro_dnk) ; %These are the german market constraints
    zeros((m.num_firms_dnk+1)*m.num_pro_dnk, (m.num_firms_ger+1)*m.num_pro_ger), kron(eye(m.num_pro_dnk), J_m_dnk) ]; 

%Fixed-Effects pattarn for Denmark - First num_firms_dnk firms affect DNK
%constraints, none impact the adding up constraint....
FE_DNK = [ [eye(m.num_firms_dnk) zeros(m.num_firms_dnk, 1) ]; zeros(m.num_firms_ger - m.num_firms_dnk, m.num_firms_dnk+1)];

%For GERMANY, all firms are active in german projects, zeros for adding up
%constraints
FE_GER = [eye(m.num_firms_ger) zeros(m.num_firms_ger,1)];

%Rows of Jacobian pertaining to fixed effects...
J_FE = [repmat(FE_GER, 1, m.num_pro_ger) repmat(FE_DNK, 1, m.num_pro_dnk)];

%Rows of Jacobian pertaining to Beta1 (dist) effect every constraint but 
%the adding up constraints...Beta2 (country) only effects foreing firms...


J_BETA_GER = [ones(1,m.num_firms_ger) zeros(1,1);...  This covers betaLogDist
              ones(1,m.num_for_ger) zeros(1,m.num_firms_ger-m.num_for_ger+1)]; %...  This covers betaFor
%              zeros(1,5) ones(1,5) 0]; %This covers betaState (assume all out of state) (last zero is adding up). 
              
          
J_BETA_DNK = [ones(1,m.num_firms_dnk) zeros(1,1);...  This covers betaLogDist
              zeros(1,m.num_firms_dnk-m.num_for_dnk) 1 0];% This covers betaFor
%              zeros(1,7)];... %German State dummy parameter (WE KNOW THAT state_mat_dnk == zeros) 

J_BETA = [repmat(J_BETA_GER,1, m.num_pro_ger), repmat(J_BETA_DNK, 1, m.num_pro_dnk)];

Jac_Pattern = [JJ; J_FE; J_BETA]';

%% Try 10 different set of starting values

r0 = [rho_ger;rho_dnk];

%Storage for the multi-start
RUNS = 3;
starts = zeros(size(x_0,1),RUNS);
sols = zeros(size(x_0,1),RUNS);
logLike = zeros(RUNS,1);
flag = -ones(RUNS,1);
out_status = cell(RUNS,1);


for i=1:RUNS

%% Solving Model for a given set of parameters...
tic;

beta0(1:m.num_firms_ger,1) = -4 + 4* rand(10,1);  % ~ U(-4,0)
beta0(m.num_firms_ger+1:m.num_firms_ger+2,1) = [ .1; 1 ] + rand(2,1);%[-.5+ rand(1,1); -2 + 2.5* rand(1,1)]; %% distance ~U(-.5,.5) , border~ U(-2,.5) 
rho = r0;
% Jac_Pattern_rho = JJ';
% LB = [zeros(size(rho_ger,1)+size(rho_dnk,1),1)];
% UB = [ones(size(rho_ger,1)+size(rho_dnk,1),1)];
% ktropts = optimset('DerivativeCheck','on','Display','iter',...
%           'GradConstr','on','GradObj','on','TolCon',1E-6,'TolFun',1E-6,'TolX',1E-6,'JacobPattern',Jac_Pattern_rho);
% [rho, FVAL, EXITFLAG, OUTPUT] = ktrlink(@(x_0) dummy_objective(x_0), r0, [],[],[],[],LB,UB,@(x_0) model_constraints(x_0, beta0,m),ktropts,'knitro.opt');
% %[rho, FVAL, EXITFLAG, OUTPUT] = fmincon(@(x_0) dummy_objective(x_0), r0, [],[],[],[],LB,UB,@(x_0) model_constraints(x_0, beta0,m),ktropts);



      
      
%% Optimization
x_0 = [rho;beta0];
starts(:, i) = x_0;

ktropts = optimset('DerivativeCheck','on','Display','iter',...
          'GradConstr','on','GradObj','on','TolCon',1E-6,'TolFun',1E-6,'TolX',1E-6,'JacobPattern',Jac_Pattern);

LB = [zeros(size(rho_ger,1)+size(rho_dnk,1),1);ones(length(beta0),1)*(-1000)];
UB = [ones(size(rho_ger,1)+size(rho_dnk,1),1);ones(length(beta0),1)*(1000)];      

[betahat,FVAL,EXITFLAG,OUTPUT] = knitromatlab(@(x_0) likelihood_func(x_0,m),x_0,[],[],[],[],LB,UB,@(x_0) eqm_constraints(x_0,m),[],ktropts,'knitro.opt');

disp(sprintf('Completed iteration %d\n', i));
toc;
sols(:,i) = betahat;
logLike(i) = FVAL;
flag(i) = EXITFLAG;
out_status{i} = OUTPUT;
save('mstart_output', 'sols', 'starts', 'logLike', 'flag', 'out_status');
% to check analytical gradients: Note that ktrlink does not provide a
% proper derviative check!
%[betahat,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x_0) likelihood_func(x_0,m),x_0,[],[],[],[],LB,UB,@(x_0) eqm_constraints(x_0,m),ktropts)

end

%%
%Here we save only the true point estimates in their own file for use
%later...in practice each multi-start run has converged to the same place,
%but just in case...
[ll, idx] = min(logLike);
betahat = sols(:,idx);
%Here we should check that this point is actually a locally optimal
%solution...
if (flag(idx) ~= 0) 
    disp(sprintf('WARNING! Exit Flag indicates this point is not locally optimal.  Flag = %d\n',flag(i)));
end
thetahat = betahat(m.num_pro_ger * (m.num_firms_ger + 1)+m.num_pro_dnk * (m.num_firms_dnk + 1)+1:end,1);
rhohat = betahat(1:m.num_pro_ger * (m.num_firms_ger + 1)+m.num_pro_dnk * (m.num_firms_dnk + 1));

save('est_point', 'betahat', 'thetahat', 'rhohat', 'll');

%% Standard errors

% 1. for k=1:length(thetahat)
%    r(thetahat(k)) = soles model thetahat+h*eps(k)
% 


%Step size for numerical derivatives..make sure step is more than contraint
%tolerance!
h = 1e-5; 
%Split entry probabilitites and paramerters out of betahat...


%Vector of which firms "win" projects (ordered in the same way as rho)
Y = [ reshape(m.y_matrix_ger', m.num_pro_ger * (m.num_firms_ger + 1), 1) ; ...
      reshape(m.y_matrix_dnk', m.num_pro_dnk * (m.num_firms_dnk + 1), 1) ];
Y = logical(Y);
  
%Number of projects (and number of winners)
N = m.num_pro_ger + m.num_pro_dnk;
K = length(thetahat);

%Allocate space for vector of partials...
dr_dth = zeros(length(rhohat), length(thetahat)); %The matrix of partials of rho with respect to thetahat

      
%Take numerical derivatives of all the rhos with respect to the parameters
%(thetahat) 
for k = 1: K
    th = thetahat;
    th(k) = th(k) + h;
    
    Jac_Pattern_rho = JJ';
    LB = zeros(size(rhohat));
    UB = ones(size(rhohat));
    ktropts = optimset('DerivativeCheck','on','Display','iter',...
              'GradConstr','on','GradObj','on','TolCon',1E-6,'TolFun',1E-9,'TolX',1E-6,'JacobPattern',Jac_Pattern_rho);
    [rk, FVAL, EXITFLAG, OUTPUT] = knitromatlab(@(x_0) dummy_objective(x_0), rhohat, [],[],[],[],LB,UB,@(x_0) model_constraints(x_0, th,m),[],ktropts,'knitro.opt');
    dr_dth(:,k) = (rhohat - rk)./h;
    if (EXITFLAG ~= 0) 
    disp(sprintf('WARNING! Model did not solve in Numerical Derivative Loop!!!  Flag = %d\n',flag(i)));
end
end
%    
%Calculate the score for each project. 
%Only the rhos of the winning firms enter the likelihood...Y is a
%logical vector that records which firm won which project.
score = dr_dth(Y,:) ./ repmat(rhohat(Y), 1, K);

%This loop creates a three diminsional matrix where each NxKxK where
%each n in the first dimension indexes the outer product of row n of
%score. 
op = zeros(N, K, K);
for i = 1:N
    op(i, :, :) = score(i,:)'*score(i,:); %The outer product...
end

% Outer product approximation of the Variance-Covariance Matrix,
% This is (AVar/N) following Hayashi's notation...
Sig_hat = inv(squeeze(sum(op, 1)));
se = diag(Sig_hat).^(.5);
    
save('est_point','Sig_hat', 'se', '-append');    
    
 diary off   