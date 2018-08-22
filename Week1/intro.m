%% Class 1: Introduction to the course and MATLAB


%% Big Picture of the Course
%
% * Learn and develop computational skills that will help you with your research.
% * Knowlwdge is applicable across fields.
% * Traditional structural fields like IO, Trade, and Macro. 
% * Theory -- simulate calibrated models and graph results.
% * Metrics -- Monte Carlo simulations. Code own estimator. 
% * Even you you don't use computational methods too much, you will have to referee papers that do. 
%

%% Why do we need computational methods?
%
% An Example:
% Consider a demand fucntion |q=p^{-0.2}|. We could easily compute the 
% quantity for any given price with a handheld calculator. We could also 
% solve for price given a quantuty.
%
% Now consider a slightly different demand: |q = 0.5p^{-0.2} + 0.5p^{-0.5}|. 
% We could easily compute demand for a given price, but if you wanted to 
% know what price cleared the market for |q=2|, no closed form inverse 
% demand function exists. We need computational methods to find the 
% price that clears the market for a given quantity. 


%% MATLAB
% MATLAB is a high-level language for numerical analysis. Code can be written
% concisely, and basic code can be easily readable by anyone who knows basic
% computer programming and linear algebra. We will go over some basic uses
% of MATLAB, but this is not a comprehensive overview. See the appendix of
% MF for a better overview, but the best way to learn is to do things 
% yourself on the homework. 
%
% MATLAB is one of many languages used to do numerical and statistical work
% in economics. Other popular choices are R, python (with numpy), and
% sometimes julia or C or even (shudder) FORTRAN. These all have pros and cons. 
% We're using MATLAB because it is user friendly and something of an industry 
% standard. It also is very well documented relative to open source
% options.
% Most importantly, all our sample code is written in MATLAB. You are free to
% program in whatever you want, just don't ask me to debug it. 
%

%% Assignment

a = 3

%%
% Terminate lines with semicolon to suppress output

b = 2;

%% 
% Matrix assignment

A = [1 2 3; 4 5 6 ; 7 8 9]

%%
% Or

B = [1 2 3;
	4 5 6;
	7 8 9]

%%
% Quick Vector of an index

c = 1:10

%%
% is the same as

d = 1:1:10

%%
% is the same as

e = linspace(1,10,10)

%%
% Assign vectors to matrices

C = [c;d;e]



%% Basic Operations 

b = [1; 2; 3]

a+b

z = a+b

A = magic(3)

D = A+B

%%
% Matrix Multiplication
A*B

%%
% Array multiplication (element-wise)
A.*B

%%
% Transpose 
A'.*B

%%
% Inverse 

inv(A)


%% Built in functions
%
% There are many built in functions. No need to add libraries or packages. Just
% google what you need and there is likely a package already in MATLAB. Some
% built in functions work well for only specific applications, and sometimes
% there are much better thrid-party packages available.
%
% We will make use of an added library called COMPECON Toolbox, <http://www4.ncsu.edu/~pfackler/compecon/toolbox.html>. 
% This should 
% already be loaded on the lab computers in the building. It is open source, so  
% download for your own personal computer. Put the folder in somewhere and 
% tell MATLAB to look for the packages by adding |addpath('path of
% folder')|
% to the top of your script or function
%
% Useful buil-in functions: |exp|, |inv|, |diag|, |ln|, |abs|, |size|, |rand|, |zeros|
% etc.  
%
%
% How to find a function? Use the documentation. Always available on the
% menu bar or you can just google what you want with `MATLAB' in the search string. 

%% Indexing

A

%%
% Refer to the first element of a matrix:
A(1)

%%
% Refer to a specific element (row X column):
A(2,3)

%%
% Refer to a vector in a matrix:
A(1:3,1)

%%
% Refer to a matrix in a matrix:
A(1:3,:)

%% Conditional statements
%
% You can execute blocks of codes conditionally
%

a=0;
if a<3
    disp('a is less than 3')
else
    disp('a is not less than 3')
end


%% Scripts and Functions
%
% A script is a file with a set of commands. If you execute the script it will
% execute those commands in order. There is no input, and the function does not
% return any output, although it can create output by printing things to screen
% or saving files. 
%
% A function is a file that is defined to take input(s) and give output(s). We 
% will often define user functions. 
%
% For example, let's say I was working on a problem that involved a demand 
% function $q = e^{-p}$. I can write a function in a separate file called
% |demand.m| for this:

%%
%
% <include>demand.m</include>
% 

%%
% Then we can call |demand.m| to evaluate any given |p|.

demand(3)

%%
% We can assign the output of a function to a variable in the local space 
% (for example in a script) 
q1 = demand(3)

p = 1:1:10;
Q = demand(p)


%% Loops
Qloop = zeros(1,10);
Qlen = length(Qloop);

for ix = 1:Qlen
	Qloop(ix)=demand(ix);
end

disp('This is the vector we just filled')
disp(Qloop)

%% No Types
% Unlike lower level languages, no need to preallocate types. Also, you can 
% write over variables, even with different types. Before we defined a variable 
% |c|. We can just write over it (this is not neccesarily a good thing!).
%
% To see this, look how an innocent matrix multiplication becomes an error.
% Also, notice how MATLAB lets me gracefully handle the error using
% try-catch. 

b = [1 2 3];

out = b*C

b = 'happy'

try 
   out = b*C
catch e
    fprintf('MATLAB is not happy: %s\n', e.message);
end



%% A Note on Efficiency 
%
% MATLAB is built differently than other languages in that it is written with 
% linear algebra in mind. In many other languages you might use loops to
% construct matrices or do algebra, but in MATLAB it is best to use the matrix
% and array operators. Also, MATLAB accesses memory differently than lower 
% level langauages like C/Fortran. Generally, you want to avoid loops. 
% Note: There is a recently developed language called *Julia* that has the 
% user-friendly features and syntax of MATLAB, but performs similarly to 
% C/Fortran. 


%% Solving Linear Equations
%
% There are many methods for solving a system of linear equations, 
% |Ax=b|. In general, MATLAB understands the properties of your system of
% equations
% and uses a method that is most efficient, yet still accurate. In general, 
% MATLAB will use an L-U factorizaiton algorithm, but will also adjust the
% matrix if it is ill-conditioned. 
%
% The most efficient way to solve a system of linear equations:

%%
% Solve system of linear equations
b = [1; 2; 3];
x = A\b

%%
% Above, MATLAB knows what is going on with the equation, and will perform 
% LU factorization (or another method if applicable) to solve for x. This is not
% the same as 

y = inv(A)*b;

%%
% although they will likely result in the same answer for a well conditioned 
% |A| matrix.

%%
% Always beware when inverting matrices. Large matricies with a wide range in 
% scaling, or with columns that are nearly correlated, can cause huge
% computational issues, and the root cause can remain hidden because MATLAB will
% technically be able to take the inverse. 



