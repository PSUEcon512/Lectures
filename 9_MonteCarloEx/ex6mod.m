% Empirical methods, Exercise 6
% based on a program by Huihui Li
% November 2011

clear
clc
delete ex6mod.out ;
diary ex6mod.out ;

A = dataset('XLSFile','robtybaer.xlsx');    % Import data
N = 650;                                    % Number of plants
T = 9;                                      % Number of years

% Reorganize dataset
pid = kron((1:N)',ones(T,1));               % Generate plant id
pid = dataset(pid);
dum = [(81:89)',eye(T)];                    % Generate year-dummy variables
dum = dataset({dum,...
    'yr','d81','d82','d83','d84','d85','d86','d87','d88','d89'});
A = horzcat(pid,A);
A = join(A,dum);

% Regressors for (a) equation (19 regressors)
Xa = {'intercept','dexp_1','dexp_2','dexp_3','d85','d86','d87','d88',...
    'd89','dcorp','dtextile','dpaper','dchemical','dmed','dother',...
    'lnage','lnw_1','lnk_1','lnpe_1'};
X.a = double(A(:,Xa));
% Regressors for (b) equation (16 regressors)
Xb = {'intercept','d82','d83','dcorp','dtextile','dpaper','dchemical',...
    'dmed','dother','lnage','lnw_1','lnk_1','lnpe_1',...
    'lnw_2','lnk_2','lnpe_2'};
X.b = double(A(:,Xb));

X.y = double(A(:,{'pid','dexp'}));      % Plant id and dependent variable
X.N = N;

% Get starting values of regression parameters from Probit models
inda = A.yr>=84;
indb = A.yr<=83;
a0 = glmfit(double(A(inda,Xa)),double(A(inda,'dexp')),'binomial',...
    'link','probit','constant','off')
b0 = glmfit(double(A(indb,Xb)),double(A(indb,'dexp')),'binomial',...
    'link','probit','constant','off')
w0 = [1.1 5 3]'            % Starting value of Omega matrix

optnm = optimset('Display','iter','MaxIter',2e5,'MaxFunEvals',5e5);
optnewton = optimset('LargeScale','off','Display','iter','TolFun',1e-5,...
                     'TolX',1e-5,'MaxFunEvals',1e7,'MaxIter',1e5);
        
%% Maximum likelihood estimation - Nelder Mead

% Optimization with 5 quadrature points

p0 = [a0;b0;w0];
[pe,fv_new,exitflag,output] = fminsearch(@(p)llk(X,p,[5 5]),p0,optnm)
b = pe;

% Optimization with 10 quadrature points

p0 = [a0;b0;w0];
[pe,fv_new,exitflag,output] = fminsearch(@(p)llk(X,p,[10 10]),b,optnm)
b = [b pe];


%  Optimization with Fminunc using bfgs - 5 quadrature points

%[pe,fv_new,exitflag,output] = fminunc(@(p)llk(X,p,[5 5]),p0,optnewton)


%% Output the results
disp(' ');
disp('Coefficient estimate of (a) equation (1984-1989)');
disp('                n = 5        n = 10 ');
for i=1:19
    fprintf('%-10s %12.4f %12.4f\n',Xa{i},b(i,1),b(i,2));
end
disp(' ');
disp('Coefficient estimate of (b) equation (1981-1983)');
disp('                n = 5        n = 10 ');
for i=20:35
    fprintf('%-10s %12.4f %12.4f\n',Xb{i-19},b(i,1),b(i,2));
end
disp(' ');
disp('Covariance matrix estimate');
disp('Quadrature point: n = 5');
disp([b(36,1) 0;b(37,1) b(38,1)]*[b(36,1) b(37,1);0 b(38,1)]);
disp('Quadrature point: n = 10');
disp([b(36,2) 0;b(37,2) b(38,2)]*[b(36,2) b(37,2);0 b(38,2)]);
disp('Log-likelihood function value');
disp('                n = 5        n = 10 ');
fprintf('           %12.4f %12.4f\n',-fv(1),-fv(2));

%% Compare to RT's results
s = repmat(b(1:19,2),1,2);
s(:,1) = sqrt(1-0.620)*s(:,1);
s(:,2) = (-7.039/s(1,2))*s(:,2);
disp(s);

diary off ;
