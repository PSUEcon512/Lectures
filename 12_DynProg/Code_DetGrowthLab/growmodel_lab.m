% ORIGINAL FILE
% r cooper 
% January 1997, modified 5/21/2001, 10/29/02
% value function iteration for
% simple non-stochastic growth model

% LAST MODIFIED by P Grieco 1/22/2014

% MODEL WITH LABOR CHOICE

clear all;% this clears variables from the workspace

% here grow.out is a diary which records aspects of the program
delete grow2.out
diary grow2.out;  % this records the output to the screen

% Section 1: Functional Forms

% the model is the simple growth model  
% preferences: u(c)=c to the power of (1-sigma) over (1-sigma)
% beta is the discount factor
% delta is the rate of capital depreciation
% the production function is: k to the power of alpha

% Section 2
% Basic Parameters

global mu alpha delta   % THESE VARIABLES ARE PASSED TO solvelab.m
                           
% global kcurrent knext %Globals are in general dangerous, especially
                        % as a means of passing non-constant data.
                        
 
       mu=0.6;   % LOG UTILITY
       beta=.9;
       alpha=.75;
       delta=0.3;

%Note that we've really dialed back on the state space size...       
N=500;

% determine the steady state capital stock 
% so here Kstar is something we solved for analytically

%This one is now an approximation since we didn't correct it for the labor
%choice
Kstar=(((1/(alpha*beta))-((1-delta)/alpha)))^(1/(alpha-1));

% put in bounds here for the capital stock space
% you may want to adjust these bounds 
Klo=Kstar*.1;     % Note that with endogenous labor steady state value of k will be
Khi=Kstar*0.5;    % LOWER 
step=(Khi-Klo)/N;

% now construct the state space using the step starting with Klo and
% ending with Khi.  by construction, the state space includes Kstar

% now create the state space as a row vector
% look up this stuff in matlab so you know how to create vectors, matrices etc.

K=Klo:step:Khi;

% n is the true length of the state space

n=length(K);

colones=ones(n,1);
rowones= colones';
I = K'*rowones;

% Section 4: SETUP
% this section has two parts. first we obtain an initial guess of the value
% function from the single period problem. 
% Second we start the value function iteration routine.
% some initial matrices

% From the production function, Kalpha is the vector of output levels

V=ones(1,n); % initial guess

% Iterations: using this first guess at V iterate 

T=200;

% FIRST WE CALCULATE U for all possible combinations
% of k and k'.
% This is done BEFORE the value function iterations
%% Solve for Labor Choice
disp('STARTING SOLUTION TO LABOR EQUATIONS');
tic

    U=ones(n,n)*(-10000);  % NOTE that we set U=-10000
    labd=zeros(n,n);
    guess=0.5; 
    
     for i=1:n;
        for j=1:n;
  
          kcurrent=K(i);   % we have to find the optimal labor supply decision for all (k,k')
          knext=K(j);
          
             
          if  ( ( ( kcurrent^alpha)*(0.9^(1-alpha)) ) - knext  +  (1-delta)*kcurrent ) >0
                 
          % NOTE that when you work 0.9 of your time
          % you do not have positive consumption
          % we do not solve for labor decision
          % then U is left at -10000 and we will not pick that k' option
          % this is a cheap way to avoid negative utility which
          % is a problem with log utility
          
          lab=fzero(@(x) solvelab(x, kcurrent, knext),guess);   % solvelab.m is a function file that defined the FOC for labor
                                      % fzero is a build in MATLAB function
                                      % note that all global variables are passed to solvelab 
           
                              
                                      
          guess=lab; 
          labd(i,j)=lab;              % optimal decision
           
          C =  ( ( ( kcurrent^alpha)*(labd(i,j)^(1-alpha)) ) - knext  +  (1-delta)*kcurrent );
            
          % consumption
             
          UC=log(C);            
          UL=log(1-labd(i,j));
          U(i,j) = (1-mu)*UC + mu*UL; % utility

          end
      end
      end
        
      U=U';
      labd=labd';

  
  %plot(labd(:,1))
  %pause
  %close
  %plot(labd(1,:))
  %pause
  %close
  
  toc
%% Value Function Iteration

%here we set our tolerance for convergence; 
%play with this to see how the results depend on tolerance

toler=.001;

% now keep doing the code from above as we loop and loop and..

% here T is the maximal number of iterations

disp('STARTING VALUE FUNCTION ITERATION');
tic

for j=1:T;
    
% use the last iteration to create a matrix of future values
% where the future capital stock is a row

% work carefully through these matrices!


% use the guess on V here to construct the matrix of payoffs

w=ones(n,1)*V;
q=w';

% this is where q which came from w which came from V comes back in

r = U + beta*q;

% CHECK HOW r looks like
% disp('iteration')
% j

v=max(r);

%plot(v)
%pause

%  now calculate the difference between the last (V) and current (v)
% guesses of the value functions

diff = sum(abs(V-v));
diffnorm = norm(V-v);
result=[j,diff,diffnorm] ;
disp(result)

% test for diff significant difference

if abs(diff) <= toler
   break
else

% so here if the diff exceeds the tolerance, then do a replacement
% and go back to the top of the loop
     V=v;
 end
end
toc

%% DIAGNOSTICS
%CHECK how v looks like
% plot(K,V)


% now that the program has converged.  You will have 
% the vector of values explicity: this is v.
% in addition, you want to know the policy vector
% use the max operator (check the Matlab book) to figure out which row was chosen for each
% level of the capital stock

[R,m]=max(r);

% now build a vector for future capital using the m-vector.

Kprime=[];
labdec=zeros(n,1);

% now loop through the m-vector using the I matrix

for i=1:n
inv=I(m(i),i);
Kprime=[Kprime inv];
labdec(i)=labd(m(i),i);
end

% here we plot the policy function against the 45 degree line

figure(1)
plot(K,Kprime,'r-');
hold on;
plot(K,K,'k:');% so now we have a nice 45 degree line
xlabel('current capital')
ylabel('capital')
legend('policy function','current capital',0)

% use this plot to look for steady state and investment patterns
figure(2)
plot(K,Kprime-K)
xlabel('current capital')
ylabel('net investment')
% this command will plot the index of the policy function
%plot(m); % 

%We also want to see the labor policy function
figure(3)
plot(K, labdec);
xlabel('current capital')
ylabel('labor')

% simulation of transition dynamics

   P=100; % arbitrary length for transition path
   capind=ones(1,P);
   captran=ones(1,P);
   capind(1)=(3); % arbitrary starting point in terms of the index
   captran(1)= K(capind(1));
   
   for t=2:P
   capind(t)=(m(capind(t-1)));% so follow evolution in index space   
   captran(t)=K(capind(t)); % follow evoluation in capital space
	end
figure(4)    
plot(captran)
xlabel('time period')
ylabel('capital stock')

  disp('the program has now ended')
