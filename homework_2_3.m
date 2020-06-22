close all
clc

load('reaction_diffusion_big.mat');
%%

[m,n,k]=size(u); % x vs y vs time data
% 
% % % Plots the video
% % for j=1:k
% % A=BZ_tensor(:,:,j);
% % pcolor(A), shading interp, pause(0.01)
% % end
% 
% 
tempu = zeros(m*n,k);
tempv = zeros(m*n,k);
for jj=1:k
    tempu = u(:,:,jj);
    uc(:,jj) = tempu(:);
    tempv = v(:,:,jj);
    vc(:,jj) = tempv(:);
end
x=[uc;vc];
%
% 
%% SVD OF THE DATA
[Ux,Sx,Vx]=svd(x,'econ');

plot(diag(Sx)/sum(diag(Sx)),'o')
% 
%% DMD with projection in the low rank subspace of the first 100 snapshots
t_train=100;
x1 = x(:,2:t_train);
x2 = x(:,1:t_train-1);
% 
rank = 10;
dt = 1;
% 
% Step 1
[U2,Sigma2,V2] = svd(x2,'econ');
U = U2(:,1:rank);
V = V2(:,1:rank);
Sigma = Sigma2(1:rank,1:rank);
% 
% Step 2
A_tilde = U'*x1*V/Sigma;
% 
% Step 3
%%
[W,D] = eig(A_tilde);

figure()
plot(real(diag(D)),imag(diag(D)),'o')
sgtitle(strcat('Eigenvalues',{' '},'Of',{' '},'A-Tilde'))
xlabel('Real'),ylabel('Imaginary')
% % 
% % Step 4
% Phi = x1*V/Sigma*W;
% lambda = diag(D);
% omega = log(lambda)/dt;
% y0 = Phi\x(:,1);
% % 
% plot(real(diag(Sigma))./sum(real(diag(Sigma))),'o')
%% FUTURE PREDICTIONS IN THE LOW RANK SPACE 
%prediction time k_future
k_future=k;
x_future=zeros(size(x));
x_future(:,1)=x(:,t_train);
x_tilde=U'*x;
x_tilde_future=zeros(size(x_tilde));
x_tilde_future(:,1)=x_tilde(:,t_train);
for i=2:k_future
    x_tilde_future(:,i)=A_tilde*x_tilde_future(:,i-1);
end
% back projection in the original data space
x_future=U*x_tilde_future;
% %% PLOT TO COMPARE
% t2=t+5;
% figure()
% plot(t,x(1,:),'-r', t2, x_future(1,:),'-k')
% figure()
% plot(t,x(131077,:),'-r', t2, x_future(131077,:),'-k')
% figure()
% plot(t,x(262144,:),'-r', t2, x_future(262144,:),'-k')
% figure()
% plot(t,x(262145,:),'-r', t2, x_future(262145,:),'-k')
% figure()
% plot(t,x(393221,:),'-r', t2, x_future(393221,:),'-k')
% figure()
% plot(t,x(524288,:),'-r', t2, x_future(524288,:),'-k')
%% BACK TO SNAPSHOT FORMAT
xu=x_future((1:262144),:);
xv=x_future((262145:524288),:);
u_future=zeros(m,n,k_future);
v_future=zeros(m,n,k_future);
for ii=1:k_future
    tempu=reshape(xu(:,ii),[m,n]);
    u_future(:,:,ii)=tempu;
    tempv=reshape(xv(:,ii),[m,n]);
    v_future(:,:,ii)=tempv;
end
%% video


% Plots the videos
for j=1:k
A=u(:,:,j);
pcolor(A), shading interp, pause(0.01)
end

for j=1:k_future
B=u_future(:,:,j);
pcolor(B), shading interp, pause(0.01)
end
%% COMPARISON BETWEEN ACTUAL u DATA AND LOW RANK DMD
figure() 
% snapshot 100
subplot(3,2,1)
A=u(:,:,101);
pcolor(A), shading interp
subplot(3,2,2)
B=u_future(:,:,1);
pcolor(B), shading interp

% snapshot 150
subplot(3,2,3)
A=u(:,:,151);
pcolor(A), shading interp
subplot(3,2,4)
B=u_future(:,:,51);
pcolor(B), shading interp

% snapshot 200
subplot(3,2,5)
A=u(:,:,201);
pcolor(A), shading interp
subplot(3,2,6)
B=u_future(:,:,101);
pcolor(B), shading interp

%% COMPARISON BETWEEN ACTUAL v DATA AND LOW RANK DMD
figure() 
% snapshot 100
subplot(3,2,1)
A=v(:,:,101);
pcolor(A), shading interp
subplot(3,2,2)
B=v_future(:,:,1);
pcolor(B), shading interp

% snapshot 150
subplot(3,2,3)
A=v(:,:,151);
pcolor(A), shading interp
subplot(3,2,4)
B=v_future(:,:,51);
pcolor(B), shading interp

% snapshot 200
subplot(3,2,5)
A=v(:,:,201);
pcolor(A), shading interp
subplot(3,2,6)
B=v_future(:,:,101);
pcolor(B), shading interp

%% PLOT u(t) and v(t) TO COMPARE CENTER point
%xpos=[1 64 128 192 256 256 256 256 256];
%ypos=[1 64 128 192 256 192 128 64 1];
xpos=[1 128 256];
ypos=[1 128 256];
t2=t+5;
figure()
for iii=1:3
    u_point=u(xpos(iii),ypos(iii),:);
    u_point=u_point(:);
    v_point=v(xpos(iii),ypos(iii),:);
    v_point=v_point(:);
    u_future_point=u_future(xpos(iii),ypos(iii),:);
    u_future_point=u_future_point(:);
    v_future_point=v_future(xpos(iii),ypos(iii),:);
    v_future_point=v_future_point(:);
    subplot(1,3,iii)
    plot(t,u_point,'-r', t2, u_future_point,'--r', t,v_point,'-b', t2, v_future_point,'--b')
    xlabel('Time')
    legend('u','u_D_M_D','v','v_D_M_D')
    title(strcat('x=',num2str(xpos(iii)),{' '},'y=',num2str(ypos(iii))))
end
    

% figure()
% plot(t,x(131077,:),'-r', t2, x_future(131077,:),'-k')
% figure()
% plot(t,x(262144,:),'-r', t2, x_future(262144,:),'-k')
% figure()
% plot(t,x(262145,:),'-r', t2, x_future(262145,:),'-k')
% figure()
% plot(t,x(393221,:),'-r', t2, x_future(393221,:),'-k')
% figure()
% plot(t,x(524288,:),'-r', t2, x_future(524288,:),'-k')