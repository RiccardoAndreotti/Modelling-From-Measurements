clear 
close all
clc

load('lynxhare.mat')

plot(years, hare,'-or', years, lynx,'-ob')
ylabel('Population [Thousands]') 
xlabel('Year') 
legend('Hare','Lynx')
title('Canadian lynx and snowshoe hare populations')

% Create data matrix
X = [hare'; lynx'];

x = hare';
y = lynx';

L=length(x);

% Calculate central finite differences
xdot = zeros (1,L-2);
ydot = zeros (1,L-2);

for i = 2:L-1
    xdot(i-1) = (x(i+1)-x(i-1))/2;
    ydot(i-1) = (y(i+1)-y(i-1))/2;
end

% Shorten the data to optimize Lokta Volterra parameters
xs = zeros (1,L-2);
ys = zeros (1,L-2);

for i = 1:L-2
    xs(i) = x(i+1);
    ys(i) = y(i+1);
end

%ODE REGRESSION

% Simple library of functions
% A=[x1s x2s x1s.^2 x1s.*x2s x2s.^2 x1s.^3 (x2s.^2).*x1s x2s.^3];

% This library contains the required terms
A=[xs' ys' xs.^2' (xs.*ys)' (ys.^2)' (xs.^3)' ((xs.^2).*ys)' ((ys.^2).*xs)' ys.^3'];

% Try the following 3 methods for computing xi1 and xi2:

%csix=A\xdot';
%csiy=A\ydot' ;
%csix=pinv(A)*xdot.';
%csiy=pinv(A)*ydot.';
csix=lasso(A,xdot,'Lambda',0.5);   % lambda determines how much sparsity you want
csiy=lasso(A,ydot,'Lambda',0.5);


figure(3)
title('Loadings')
subplot(2,1,1), bar(csix)
title('\xix')
subplot(2,1,2), bar(csiy)
title('\xiy')
%% 

dt=1.0;
t=1:dt:30; 
x0=[20; 32]; % initial conditions


[t,y]=ode45('sparse_mod',t,x0);

figure()
plot(years, y(:,1),'-or', years, y(:,2),'-ob')
ylabel('Population [Thousands]') 
xlabel('Year') 
legend('Hare','Lynx')
title('Lynx and hare populations according to the model')
