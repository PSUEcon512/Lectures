% ORIGINAL FILE
% r cooper 
% January 1997, modified 5/21/2001, 10/29/02
% value function iteration for
% simple non-stochastic growth model

% MODIFIED by N Guner 1/17/2005
% MODIFIED AGAIN by P Grieco 11/24/2010

clear all;% this clears variables from the workspace

% here grow.out is a diary which records aspects of the program
delete grow.out
diary grow.out;  % this records the output to the screen

% Section 1: Functional Forms

% the model is the simple growth model  
% preferences: u(c)=c to the power of (1-sigma) over (1-sigma)
% beta is the discount factor
% delta is the rate of capital depreciation
% the production function is: k to the power of alpha

% Section 2
% Basic Parameters

     sigma=1;
     beta=.9;
     alpha=.75;
     delta=0.3;

%%
% Section 3: Spaces
% here you need to specify the space that the capital stock lives in
% you must make an educated guess here so that the capital stock space
% includes the steady state
% input the size of the state space
% this will be the number of points in the grid

N=1000;

% determine the steady state capital stock 
% so here Kstar is something we solved for analytically

Kstar=(((1/(alpha*beta))-((1-delta)/alpha)))^(1/(alpha-1))

% put in bounds here for the capital stock space
% you may want to adjust these bounds 

Klo=Kstar*.9;
Khi=Kstar*1.1;
step=(Khi-Klo)/N;

% now construct the state space using the step starting with Klo and
% ending with Khi.  by construction, the state space includes Kstar

% now create the state space as a row vector
% look up this stuff in matlab so you know how to create vectors, matrices etc.

K=Klo:step:Khi;


% n is the true length of the state space

n=length(K);

%%
% Section 4: PREP FOR VALUE FUNCTION ITERATION
% this section has two parts. first we obtain an initial guess of the value
% function from the single period problem. 
% Second we start the value function iteration routine.
% some initial matrices

% From the production function, Kalpha is the vector of output levels

Kalpha=K.^alpha;
disp(['The dimensions of Kalpha are ' num2str(size(Kalpha)) ])
% ytot is then total output available

ytot=Kalpha+(1-delta)*K; 

% create a column of n ones 

colones=ones(n,1);

% create an nxn matrix here
% Could also use repmat...
s = colones*ytot;

% now take the transpose

S=s';

% each column of S is ytot
% we created all of this so we can work with matrices and avoid loops

% first guess of the value function

% here describe the level of utility associated with each element of the
% matrix S.  we use this below as our first guess for the dynamic programming
% problem.

% that is this is our initial guess on V

if sigma==1
V=log(S);
else
V=S.^(1-sigma)/(1-sigma);
end
disp(['The dimensions of initial V are ' num2str(size(V)) ])
Vcol1=V(:,1); 

% calculate investment here; I plays the role of the future capital stock
% this isn't great notation so be careful here....

rowones= colones';
% so I here is matrix where each column is K
I = K'*rowones;
% here J is a matrix where each row is K
J= colones*K;

% consumption: the flow of output plus undepreciated capital less investment
% you may want to work this out to see that it is correct using the matrices

C = (J.^alpha)-I +(1-delta)*J ;
%ct=(J.^alpha) +(1-delta)*J ;
%ctcol1=ct(:,1);
%Icol1=I(:,1);
%Ccol1=C(:,1);
%Ccollast=C(:,n);

% current utility, current state as cols, future capital as rows
if sigma==1
U=log(C);
else
U=(C.^(1-sigma))/(1-sigma);
end


% create the first iterate of the payoffs to the dynamic problem
% using the U just constructed and V as described above
% If I am going from state "J" to state "I" I get U utility today
% then my continutation utility (guess) is V, which depends only on
% Where I go, (state "I"). 
r = U + beta*V ;

% Here's what optimization over a discrete set of choices (which state to 
% go to) looks like: it is a column max so that we are effectively choosing next periods
% capital for each value of today's capital

V=max(r);
%Now V is just a vector, we'll need to replicate it to carry out the 
%next iteration...
disp(['The dimensions of V are ' num2str(size(V)) ])
Vrow1=V(1,:); 


% Iterations: using this first guess at V iterate 

T=300;

%here we set our tolerance for convergence; 
%play with this to see how the results depend on tolerance

toler=.001;

% now keep doing the code from above as we loop and loop and..

% here T is the maximal number of iterations

%%
%%SECTION 5: VALUE FUNCTION ITERATION
for j=1:T;
% use the last iteration to create a matrix of future values
% where the future capital stock is a row

% work carefully through these matrices!


% use the guess on V here to construct the matrix of continuation payoffs
q = repmat(V',1,n);
%Anothe way of doing the same thing...
%w=ones(n,1)*V;
%q=w';

% R are the choice-specific continuation values given V...

r = U + beta*q ;

% CHECK HOW r looks like

% subplot(2,1,2)
% plot(K,r(:,1))
% subplot(2,1,2)
% plot(K,r(:,50))

%The optimization to get a new candidate for v
v=max(r);

%  now calculate the difference between the last (V) and current (v)
% guesses of the value functions

diff=sum((abs(V-v))) ;
diffnorm=norm(V-v);
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

disp('  convergence achieved  ')


%%
%SECTION 6: Output...

%CHECK how v looks like
figure(1) 
plot(K,V)
xlabel('current capital')
ylabel('Value Function')
%disp('  Current capital ')
%disp(K) 
%disp('  Value function  ')
%disp(V)

% now that the program has converged.  You will have 
% the vector of values explicity: this is v.
% in addition, you want to know the policy vector
% use the max operator (check the Matlab book) to figure out which row was chosen for each
% level of the capital stock

[R,m]=max(r); 

% now build a vector for future capital using the m-vector.
Kprime=K(m);



figure(2)
plot(K,Kprime,'r');
hold on;
plot(K,K);% so now we have a nice 45 degree line
xlabel('current capital')
ylabel('next period capital')
%legend('policy function','current capital',0)

% use this plot to look for steady state and investment patterns
figure(3)
plot(K,Kprime-K)
xlabel('current capital')
ylabel('net investment')
% this command will plot the index of the policy function
%plot(m); % 

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
%disp ('captran')
%disp (captran) 
%disp ('capind  ')
%disp (capind)

  disp('the program has now ended')
