%% Gravity model Estimation
% Code written by Felix Tintelnot, March 2010
% Modifications by M. Roberts, Sept 2010; Markdown by C Murry, 2015

%% Outline of code
% 
% * Guess a P vector (recall P is the price index for each country, which is a solution of a system of non-linear equations)
% * Given P vector estimate the gravity model
% * Using trade cost estimates, solve for new P vector using root-finding
% technique
% * Repeat using new P vector


%% Top Matter
clear;
delete main-est.out 
diary main-est.out
addpath('../CEtools/');

warning('off','all')
warning

disp(' ');
disp('Standard A&vW gravity estimation ');
disp(' ');



%% Load Data and Set Parameters
% Data is in following order: 

%%
% * dist exp_gdp imp_gdp trade_value language
% * history EU NAFTA ASEAN
% * Missing values are denoted by 99999
% * The bilateral trade flows, xij, and distances, dij, where
% * i=exporter, j=importer are ordered (11,12,...,1N,21,22....,2N,...NN)'

load trdata.mat;

trade   = trdata(:,4);    %Bilateral trade flows (cif values)
dis     = trdata(:,1);    %Distance
expgdp  = trdata(:,2);    %Exporter's GDP
impgdp  = trdata(:,3);    %Importer's GDP
lang    = trdata(:,5);    %dummy = 1 if both countries share an official language or >9% of population speak same language
hist    = trdata(:,6);    %dummy = 1 if both countryies share common history (same country in the past or colonial ties)
eu      = trdata(:,7);    %dummy = 1 if both countries are member EU27/EFTA
nafta   = trdata(:,8);    %dummy = 1 if both countries are member NAFTA
asean   = trdata(:,9);    %dummy = 1 if both countries are member ASEAN

N       = sqrt(size(trade,1));     %Number of countries
m.N     = N;
m.y     = impgdp(1:N,1);     % GDP of each country
m.sigma = 5;

%%
% Eliminate bilateral trade flows that are missing 
mis     = 0;                  
for i=1:N*N                   % Count number of missing tradeflows
    if trade(i) == 99999
mis=mis+1;
    end
end    
disp (['missing number of trade flows ' num2str(mis)])

i=1;
j=1;
id  = 1:1:N*N;        % observation number 
id=id';

nid=zeros(N*N-mis,1);  % nid contains the obs number of all non-missing obs  
while i <= N*N
    if trade(i) ~= 99999;
        nid(j) = id(i);
        j=j+1;
    end
    i=i+1;
end

%% 
% Eliminate all missing trade flows in the data.  These variables are
% used in the NL least squares regression

m.tr     = trade(nid);        
m.di     = dis(nid);
m.eg     = expgdp(nid);
m.ig     = impgdp(nid);
m.la     = lang(nid);
m.hi     = hist(nid);
m.eu     = eu(nid);
m.na     = nafta(nid);
m.as     = asean(nid);

%%
% We could also do this using indexing:
% * nid = trade~=99999


%% Start computational algorithm, standard A&vW regression

%% 
% *Initial Guesses*

P       = ones(N,1);     % Country Price index, initial value      
beta    = ones(7,1);     % Coefficients to be estimated, initial guess
Pnew    = P;
Pold    = P*5;
j       = 1;


%% 
% *Outer iteration over the Price Indicies (while loop)*

tol = 1e-11 ; 
while abs(Pnew-Pold)>tol
%disp(' ') ;    
%disp('start beta iteration number ');
%fprintf('%2i\n', j);

P       = Pnew;

%%
% * Country price index for the exporting country (Pexprel)
% * each block of N values is the price index for one exporting country
% * Country price index of the importing country (Pimprel)
% * each block on N values is the price index for the N importing countries

unitv=ones(N,1) ;
Pexprel = kron(P,unitv) ;
Pimprel = kron(unitv,P) ;
 
% drop obs with missing trade flows in the estimation
m.Pex     = Pexprel(nid);       
m.Pim     = Pimprel(nid);

%% 
% *Minimization algorithm*
%
% - compute beta that minimizes ssr each iteration uses beta from the previous iteration as start values
%
% - this is conditional on the current guess of |P|
%
% - This is a non-linear least squares problem, so use |lsqnonlin|

option1 = optimset('TolX',1e-4,'TolFun',1e-4,'MaxFunEvals',2000,'MaxIter',100,'Algorithm','levenberg-marquardt','Display','off');
[beta,resnorm,residual,exitflag,output]=lsqnonlin(@(beta) gravity_estobj(beta,m),beta,-inf,inf,option1);  

%%
% *Funtion being minimized*
% 
% <include>gravity_estobj.m</include>
%


%%
% Display Results

%disp([' finish beta iteration number ' num2str(j)]) ;
beta 
resnorm
exitflag
output

%% Calculate trade costs based on new gravity coefficients
% For a given beta vector, calculate trade costs and solve for 
% country price indexes

temp = beta(3)*(1-lang) + beta(4)*(1-hist) + beta(5)*eu + beta(6)*nafta + beta(7)*asean; 
tcprel   = (dis .^beta(2)) .* exp(temp);

%%
% domestic trade costs normalized to one
for i=1:N                  
    tcprel((i-1)*N+i)=1;
end;

% N by N matrix of trade costs
m.tc   = reshape (tcprel,N,N);

%% Solve for new implied P Index (using Broyden)
% Solving equation system for country Price indexes. 
% We use the full data to solve for country price indexes,
% and use P vector from previous iteration as start values

%%
% Set options for Broyden
optset('broyden','showiters',1000);
optset('broyden','maxit',300) ;
optset('broyden','tol',1e-8) ;

%%
% *Call Broyden*
tic
[P,pfval]  = broyden(@(P) gravity_priobj(P,m),P);  
toc
Pold=Pnew;
Pnew=P;

j=j+1;
end

%% final coefficient vector 
disp(' ');
disp('(1) Finished Standard A&vW gravity estimation');
m.sigma
betafinal=beta 
Pnew
lnp=log(Pnew)

diary off ;