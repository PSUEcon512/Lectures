


%function to approximate: f(x) = exp(sin(x*z))
%Set z = [1, 2, 20]... Controls periodicity of the function. 
z = 2;

%This generates the approximation nodes and the Phi matrix which 
%contains the chebyshev function values at the nodes
cD = 5; %Number of nodes to use..
[nodes, Phi] = setupCheb(cD,1,5);

%Calculate the function to approximate at the nodes
%fn = exp(sin(nodes)); 
fn = exp(sin(nodes*z));
%fn = exp(sin(nodes*20));

%Do the approximation recall fn = Phi*c
%Where vectors are column vectors...
c = Phi\fn;

%Now define a grid to check the approximation
%We're using column vectors for the list of points
grid = [1:.01:5]';

%Here we actually evaluate approximation at every point on the grid.
fhat = chebEval(grid, c, cD, 1, 5);

%And plot the results...
%NOTE: I include an offset so that we can actually see the two colors.
close all;
plot(grid, fhat+.01);
hold on;
scatter(nodes, fn, 50, 'filled');
plot(grid, exp(sin(grid*z)), 'r');
