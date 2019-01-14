% Stochastic Version of One Sector Growth Model
% Penn State, Econ 512

clear all;
close all;

% Parameter values

beta=0.95; % discount factor
alpha=0.3; % labor share
gamma=2;   % utility parameter
%No depreciation in this version, but it is easy to add...
  
 Z=21;              % grid number for shocks
 sigma=0.02;       % std of error term in AR(1) process
 mu=0;           % unconditional mean 
 ro=0.81;       % AR(1) parameter 
 
 

% GRID FOR CAPITAL

% step 1: using Tauchen method to find the grid points for the erogodic set
% of In(k)

N=1000;
varz=(sigma)^2;
vare=(1-ro^2)*varz;

%stde=sqrt(vare);

muk=((1-alpha)^(-1)*(log(alpha*beta)+vare/2));
vark=(ro^2)*vare/((1-ro^2)*((1-alpha)^2));
logk=muk-4*sqrt(vark):8*sqrt(vark)/(N-1):muk+4*sqrt(vark);

% step 2: approximate the expected value of capital by using the grid points
kk=exp(logk);
kkm=mean(kk);

% step 3: find the grid point for k in between 20% below and 20% above the
% expected value of k

N=100;
khigh=1.2;
klow=0.8;
kdif=khigh-klow;

k=klow*kkm:kdif*kkm/N:(khigh-(1-klow)/N)*kkm;

clear varz vare muk vark logk kk kkm; 

% GRID FOR PRODUCTIVITY SHOCKS

% grid for shock by using Tauchen's method for finite state Markov
% approximation

[prob,grid]=tauchen(Z,mu,ro,sigma);
disp(['The dimensions of prob are ' num2str(size(prob)) ])
disp(['The dimensions of grid are ' num2str(size(grid)) ])
% VALUE FUNCTION ITERATION

vinitial=zeros(N,Z);
vrevised=zeros(N,Z);
decision=zeros(N,Z);
  
invest=kron(ones(1,Z),k');
disp(['The dimensions of invest  are ' num2str(size(invest)) ])
ONE=ones(N,1);
 
%iteration

diff=1;

while diff>0.001
    
    Ev=vinitial*prob';   % find the expected value of value function
      
    for i=1:N           % for each k, find the optimal decision rule
      ci=kron(exp(grid),k(i)^alpha*ONE)-invest;  
     % disp(['The dimensions of ci are ' num2str(size(ci)) ])
      [vrevised(i,:),decision(i,:)]=max((ci.^(1-gamma))/(1-gamma)+beta*Ev);
     % disp(['In the loop vreviesed are ' num2str(size(vrevised)) ])
    end
   % disp(['The dimensions of vreviesed are ' num2str(size(vrevised)) ])
     diff=norm(vrevised-vinitial)/norm(vrevised)
     vinitial=vrevised;
    
end

% compute decision rule

derule=zeros(N,Z);

for i=1:Z
    derule(:,i)=k(decision(:,i))';
end


plot(k,derule)
hold on
plot(k,k)
xlabel('current capital');
ylabel('next period capital');
title('optimal investment decision')
%disp('PRESS ANY KEY TO CONTINUE')
%pause
%close


% TRANSITION MATRIX and LONG RUN DISTRIBUTION FOR CAPITAL

P=zeros(Z*N,Z*N);

for i=1:Z;
    for j=1:Z;
        P((i-1)*N+1:i*N,(j-1)*N+1:j*N)=prob(i,j)*(kron(ones(1,N),derule(:,i))==kron(ones(N,1),k));
        
    end
end

% Finding stationary distribution

px=ones(1,Z*N);
pinitial=px/(Z*N);

m=1;

while m>0.001
    px=pinitial*P;
m=norm(px-pinitial)/norm(pinitial);
pinitial=px;
end

probk=zeros(N,Z);

for i=1:Z;
    probk(:,i)=px((i-1)*N+1:i*N)';
end

probk=sum(probk');

figure(2);
plot(k,probk)
title('long run distribution of capital')
