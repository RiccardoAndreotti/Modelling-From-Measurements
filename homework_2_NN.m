clear all
close all
clc

maxtime = 0.6;
n_solutions = 100;
N = 100;

input = [];
output = [];
for jj = 1:n_solutions
    
    [tsave,xsave,usave,dt,dx] = KS_function(maxtime,N);
    usave = usave.';
    
    % Plot training data
    surf(xsave*dx,tsave*dt,usave.'), shading interp, title('Training data')
    xlabel('x'), ylabel('time'), zlabel('u')
    pause(0.01)
    
    % Create Input and Output data for the Neural Network
    input = [input usave(:,1:end-1)];
    output = [output usave(:,2:end)];
    
end

%%

net = feedforwardnet([200 200 200 200 200 200 ]);
net.layers{1}.transferFcn = 'tansig'; 
net.layers{2}.transferFcn = 'tansig'; 
net.layers{3}.transferFcn = 'tansig'; 
net.layers{4}.transferFcn = 'tansig'; 
net.layers{5}.transferFcn = 'tansig'; 
net.layers{6}.transferFcn = 'tansig'; 
net.trainFcn = 'trainscg';
net.trainParam.epochs=100000;
net.trainParam.max_fail = 1000;
net = train(net,input,output,'useGPU','yes');

%% test net

% Neural Net Solution #1
unn(:,1)=input(:,1);% sets the ICs as the first state in the solution

for jj=2:length(tsave) 
    unn(:,jj)=net(unn(:,jj-1));
end
    
figure(57)
subplot(1,2,2)
pcolor(xsave,tsave,unn.') ,shading interp, xlabel('x'), ylabel('time')
title('Training set 1 - NN')
subplot(1,2,1)
pcolor(xsave,tsave,[input(:,1:length(tsave)-1) output(:,length(tsave)-1)]')
shading interp, xlabel('x'), ylabel('time'),
title('Training set 1 - ODE')


% Neural Net Solution #2 - New

[tsave,xsave,usave,dt,dx] = KS_function(maxtime,N);
usave = usave.';

clear unn
unn(:,1)=usave(:,1);% sets the ICs as the first state in the solution

for jj=2:length(tsave) 
    unn(:,jj)=net(unn(:,jj-1));
end

figure(58)
subplot(1,2,1)
pcolor(xsave,tsave,usave.'), shading interp, title('New ICs - ODE')
xlabel('x'), ylabel('time')
subplot(1,2,2)
pcolor(xsave,tsave,unn.'), shading interp, xlabel('x'), ylabel('time')
title('New ICs - NN')







