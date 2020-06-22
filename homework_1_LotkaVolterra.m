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
x_dot = zeros (1,L-2);
y_dot = zeros (1,L-2);

for i = 2:L-1
    x_dot(i-1) = (x(i+1)-x(i-1))/2;
    y_dot(i-1) = (y(i+1)-y(i-1))/2;
end





% Shorten the data to optimize Lokta Volterra parameters
x_short = zeros (1,L-2);
y_short = zeros (1,L-2);
years_short = zeros (1,L-2);

for i = 1:L-2
    x_short(i) = x(i+1);
    y_short(i) = y(i+1);
    years_short(i) = years(i+1);
end


figure()
plot(years_short, x_short,'-or', years_short, y_short,'-ob', years_short, x_dot,'-og', years_short, y_dot,'-oc')
ylabel('Population [Thousands]') 
xlabel('Year') 
legend('Hare','Lynx','Hare finite derivative','Lynx finite derivative')
title('Canadian lynx and snowshoe hare populations')



% Using matrix formulation to solve the overconstrained problem and calculate the opimized parameters for x

Ax=zeros(L-2,2);

Ax(:,1)=x_short';
Ax(:,2)= -1*(x_short'.*y_short');

csi_x = pinv(Ax)*x_dot'

b=csi_x(1)
p=csi_x(2)



% Using matrix formulation to solve the overconstrained problem and calculate the opimized parameters for x

Ay=zeros(L-2,2);

Ay(:,1)= -1*(y_short');
Ay(:,2)= x_short'.*y_short';

csi_y = pinv(Ay)*y_dot'

d = csi_y(1)
r = csi_y(2)

%%
TSPAN=zeros(size(years));
for i=1:size(years)
    TSPAN(i)=i-1
end
%%
[t,yp]=ode45(@lv,TSPAN,[20;32]);
figure()
plot(years,yp)

plot(years, yp(:,1),'-or', years, yp(:,2),'-ob')
ylabel('Population [Thousands]') 
xlabel('Year') 
legend('Hare','Lynx')
title('Lynx and Hare populations according to Lotka-Volterra model')
