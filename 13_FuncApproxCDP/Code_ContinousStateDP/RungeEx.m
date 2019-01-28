%% The Runge phenomenon
% From the teaching slides of Donna Calhoun, http://math.boisestate.edu/~calhoun/

%% Introduction
% In this example, we will investigate the Runge Phenomenon and see how we might 
% be able to fix it by choosing interpolating points wisely.
% 

%% 
% As a classic example of the Runge Phenomenon, we try to interpolate the
% function f(x) = 1/(1 + 25x^2) at equally spaced points. 

clear all;
close all;

%%
f = @(x) 1./(1 + 25*x.^2);

%% 
% Plot the original function f(x) at set of equally spaced points
%

figure(1);
clf;

x = linspace(-1,1,500);
y_true = f(x);
plot(x,y_true,'r','linewidth',2);
hold on;


%%
% <html>
% Construct an interpolating polynomial of degree N at N+1 points. 
% We will use <span class="var">polyfit</span> although the problem occurs
% with Lagrange polynomial interpolation or the Barycentric form as well. 
% </html>

% Equally spaced points
N = 10;    % Degree of the polynomial we try to fit
xdata = linspace(-1,1,N+1)';
ydata = f(xdata);

p = polyfit(xdata,ydata,N);


%%
% Plot the interpolating polynomial

y_fit = polyval(p,x);

poly_10 = plot(x,y_fit,'b','linewidth',2);

plot(xdata,ydata,'k.','markersize',30);
snapnow;


%%
% We see that a single interpolating polynomial does a particularly poor
% job of interpolating near the endpoints of the interval [-1,1].  Perhaps
% we should try increasing the order of accuracy. 

N = 20;  % Increase the degree to 20 from 10. 
xdata = linspace(-1,1,N+1)';
ydata = f(xdata);

% Ignore the warning about the matrix being ill-conditioned.  If we used
% the barycentric form, we do not have these ill-conditioning issues!
p = polyfit(xdata,ydata,N);


%%
% Plot the interpolating polynomial

y_fit = polyval(p,x);

poly_20 = plot(x,y_fit,'g','linewidth',2);

plot(xdata,ydata,'k.','markersize',30);

axis([-1 1 -5 5]);
snapnow;


%%
% Let's increase one more time to see that pattern

N = 40;  % Increase the degree to 40 from 20. 
xdata = linspace(-1,1,N+1)';
ydata = f(xdata);

% plot(xdata,ydata,'k.','markersize',30);

% Ignore the warning about the matrix being ill-conditioned.  If we used
% the barycentric form, we do not have these ill-conditioning issues!
p = polyfit(xdata,ydata,N);


%%
% Plot the interpolating polynomial

y_fit = polyval(p,x);

poly_40 = plot(x,y_fit,'m','linewidth',2);

plot(xdata,ydata,'k.','markersize',30);

axis([-1 1 -10 10]);
snapnow;

%%
% But are we doing at least a better job in the middle? Yes!

axis([-1 1 0 1]);
lh = legend([poly_10, poly_20, poly_40],{'N = 10','N = 20','N = 40'});
set(lh,'fontsize',18);
snapnow;

