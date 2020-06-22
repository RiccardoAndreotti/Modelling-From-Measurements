clear all
close all
clc

load('BZ.mat');
BZ_tensor = imresize(BZ_tensor,1/8);

[m,n,k] = size(BZ_tensor); 
BZ = zeros(m*n,k);
for jj=1:k
    temp = BZ_tensor(:,:,jj);
    BZ(:,jj) = temp(:);
end

n_train = 500;
n_forecast = 50;
input = BZ(:,1:n_train-1);
output = BZ(:,2:n_train);

net = feedforwardnet([200 200 200]);
net.layers{1}.transferFcn = 'logsig'; 
net.layers{2}.transferFcn = 'logsig'; 
net.layers{3}.transferFcn = 'logsig'; 
net.trainFcn = 'trainscg';
net.trainParam.epochs=100000;
net.trainParam.max_fail = 1000;
net = train(net,input,output,'useGPU','yes');

unn = zeros(m*n,n_forecast);
unn(:,1)=output(:,end); % first frame

for jj=2:n_forecast 
    unn(:,jj)=net(unn(:,jj-1),'useGPU','yes');
end

NN_video = reshape(unn,[m n n_forecast]);

for j=1:n_forecast
    subplot(1,2,2)
    A=NN_video(:,:,j);
    pcolor(A), shading interp, axis tight,  title('NN forecast')
    subplot(1,2,1)
    A=BZ_tensor(:,:,j+n_train-1);
    pcolor(A), shading interp, axis tight,  title('Original data'),pause(0.1)
    sgtitle(strcat('frame ',{' '},num2str(j)))
end



