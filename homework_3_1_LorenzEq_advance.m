clear all
close all
clc

n_training = 1000; % Number of training trajectories
rho_training = [10 28 40];

input=[]; output=[];

for j=1:n_training  % generate training trajectories

    r = rho_training(randi(length(rho_training))); % randomly selects between the training rho values
    x0=30*(rand(3,1)-0.5);      % random ICs
    [t,y] = lorenz_function(x0,r);
    input=[input; y(1:end-1,:)];
    output=[output; y(2:end,:)];
    
    plot3(y(:,1),y(:,2),y(:,3)), hold on
    plot3(x0(1),x0(2),x0(3),'ro')
end

grid on, view(-23,18)
xlabel('x'),ylabel('y'),zlabel('z')
title('Training Trajectories - Red circle = Initial conditions')

net = feedforwardnet([10 10 10]);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas'; 
net.layers{3}.transferFcn = 'purelin'; 
net.trainFcn = 'trainscg';
net.trainParam.epochs=100000;
net.trainParam.max_fail = 1000;
net = train(net,input.',output.','useGPU','yes');



%%
load('trained_net.mat')

r = 35;
close(figure(2))
figure(2)

x0=30*(rand(3,1)-0.5); % New initial conditions

[t,y] = lorenz_function(x0,r);
plot3(y(:,1),y(:,2),y(:,3)), hold on
plot3(x0(1),x0(2),x0(3),'ro','Linewidth',[2])
grid on

% Neural Net Solution #1
ynn(1,:)=x0; % 
for jj=2:length(t) 
    y0=net(x0);  
    ynn(jj,:)=y0.'; 
    x0=y0;          
end

plot3(ynn(:,1),ynn(:,2),ynn(:,3),':','Linewidth',[2])
xlabel('x'),ylabel('y'),zlabel('z')

% --- picks middle of trajectory as IC for the NN -----
close(figure(3))
figure(3)

x0=30*(rand(3,1)-0.5); % New initial conditions

[t,y] = lorenz_function(x0,r);
plot3(y(:,1),y(:,2),y(:,3)), hold on
plot3(x0(1),x0(2),x0(3),'ro','Linewidth',[2])
grid on

% Neural Net Solution #2
ynn(1,:)=y(200,:); 
for jj=2:length(t) 
    y0=net(x0);  
    ynn(jj,:)=y0.'; 
    x0=y0;          
end

plot3(ynn(:,1),ynn(:,2),ynn(:,3),':','Linewidth',[2])
xlabel('x'),ylabel('y'),zlabel('z')


