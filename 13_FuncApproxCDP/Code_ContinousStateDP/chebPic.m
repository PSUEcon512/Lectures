
%This little script lets you visualize how orthogonal polynomials are able
%to efficiently span the entire space of continous functions from -1 o 1. 

clear all;
close all;

%Try orders 5, 10, 15, 50...
ord = 5;
X = [-1:.01:1]';
funcs = chebFuncs(X, ord, -1, 1);
plot(X, funcs);
hold on;
plot(X, zeros(size(X)))

