%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates the posterior distribution for a parameter vector
% theta supplied by the function fun3.m. The parameters are drawn 
% sequentially from the proposal distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%   construct a new data set of 500 draws from a N(2,4) distribution or load 
%   an existing data set

% mu_a  = 2 ;
% var_a = 4 ;
% cases = 500;   
% data_n24 = mvnrnd(mu_a,var_a,cases);
% save data_n24;
% f_out = log(mvnpdf(data_n24,mu_a,var_a));

 load data_n24;
 global data
 data = data_n24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First sets means and variances of the prior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



prsig      = 500*eye(2);
prmu       = zeros(1,2);
prmu(1)    = 10;
prmu(2)    = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose standard deviations for innovations in proposed theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

promu       = zeros(2,1);
prosig      =.1*ones(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Initialize the theta chain and the proposal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

propose    = zeros(1,2);
theta0 =  5*ones(1,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Construct the Markov chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T          = 50000;             % total length of the chain
iaccept    = zeros(1,2);       % counter for number of acceptances
theta(1,:) = theta0;
curr_pi    = fun3(theta0) + log(mvnpdf(theta0,prmu,prsig));

t = 2;
while t <= T;  
    j=1;
      propose   = theta(t-1,:);
    while j <=2;
      innov     = normrnd(promu(j),prosig(j));
      propose(j)= propose(j) + innov;       
      if j == 2; propose(j)= abs(propose(j)); end; % draw from half-normal for sigma
      lik       = fun3(propose);                   % log-likelihood at proposed theta
      prior     = log(mvnpdf(propose,prmu,prsig)); % log of prior at proposed theta 
      prop_pi   = lik + prior;
      delta     = prop_pi - curr_pi ;
      accprob   = min(0,delta);     
      u = unifrnd(0,1);
      u = log(u);
         if u < accprob; iaccept(j) = iaccept(j) + 1;
             theta(t,j) = propose(j);
             curr_pi = prop_pi;
         else
             theta(t,j) = theta(t-1,j);
             propose(j) = theta(t-1,j);
         end;
      j=j+1;
      end;
   t = t+1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Construct summary statistics for chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save theta ;
mtheta = mean(theta)
vtheta = var(theta)
acc_rate = iaccept./T

outfile=fopen('MCMC3.out','w');
format('compact');
fprintf(outfile,'\n');
fprintf(outfile,'SUMMARY STATISTICS\n');
fprintf(outfile,'---------------------------------------------------------\n\n');
fprintf(outfile,'  Mean of posterior mu distribution:            %f\n',mtheta(1));
fprintf(outfile,'  Mean of posterior sigma distribution:         %f\n',mtheta(2));
fprintf(outfile,'  Variance of posterior mu distribution:        %f\n',vtheta(1));
fprintf(outfile,'  Variance of posterior sigma distribution:     %f\n',vtheta(2));
fprintf(outfile,'  Acceptance rate, mu chain:                    %f\n',acc_rate(1));
fprintf(outfile,'  Acceptance rate, sigma chain:                 %f\n',acc_rate(2));
