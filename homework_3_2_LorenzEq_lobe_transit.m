clear all
close all
clc

n_training = 2000; % Number of training trajectories
n_test = 1000;
rho = 28;
n_y = 801;

input=[]; output=[];

y0 = zeros(n_training,3);

figure(1)
hold on

for jj=1:n_training  % generate training trajectories

    x0=30*(rand(3,1)-0.5);      % random ICs
    [~,y1] = lorenz_function(x0,rho);  % creates random trajectories 
    y0(jj,:) = y1(end,:);     % use end of the y1 trajectories as new ICs
    
    [t2,y2] = lorenz_function(y0(jj,:),rho);
    
    
    for k=1:n_y-1
        if y2(k,1)*y2(k+1,1)<0        
            figure(1)
            plot3(y2(k,1),y2(k,2),y2(k,3),'r*')
            input = [ input; y2(1:k,:)];
            output = [ output; flipud(t2(1:k))];
            
            break
        end
    end
end


plot3(input(:,1),input(:,2),input(:,3),'.')
grid on, view(-23,18)
xlabel('x'),ylabel('y'),zlabel('z')

net = feedforwardnet([30 30 30 30 30 30 30 30 30 30]);
net.layers{1}.transferFcn = 'tansig'; 
net.layers{2}.transferFcn = 'tansig'; 
net.layers{3}.transferFcn = 'tansig'; 
net.layers{4}.transferFcn = 'tansig';
net.layers{5}.transferFcn = 'tansig'; 
net.layers{6}.transferFcn = 'tansig';
net.layers{7}.transferFcn = 'tansig';
net.layers{8}.transferFcn = 'tansig'; 
net.layers{9}.transferFcn = 'tansig';
net.layers{10}.transferFcn = 'tansig';

net.trainFcn = 'trainscg';
net.trainParam.epochs=100000;
net.trainParam.max_fail = 1000;
net = train(net,input.',output.','useGPU','yes');


%%

%load('saved_net_time_forecast.mat')

jj=1;
test_input=[];
test_output=[];

while (jj<=n_test)  % generate training trajectories

    x0=30*(rand(3,1)-0.5);      % random ICs
    [~,y1] = lorenz_function(x0,rho);  % creates random trajectories 
    y0(jj,:) = y1(end,:);     % use end of the y1 trajectories as new ICs
    
    [t2,y2] = lorenz_function(y0(jj,:),rho);
    
    for k=1:n_y-1
        if y2(k,1)*y2(k+1,1)<0
            test_input = [ test_input; y2(k,:)];
            test_output = [ test_output; t2(k)];
            jj = jj+1;
            break
        end
    end
end

test_res=zeros(n_test,1);

for jj = 1:n_test
    test_res(jj) = net(test_input(jj,:).','useGPU','yes');
end



%% 

%load('test_res.mat')

dif =abs((test_res-test_output)./test_res);
figure()
histogram(abs(test_output))
figure()
histogram(abs(test_res-test_output),500)

%% Bar chart

toll = 0.1;
correct = zeros(n_test,1);

temp = dif (dif<toll);
correct = length(temp)






