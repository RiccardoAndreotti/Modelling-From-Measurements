clearvars -except BZ_tensor BZ
close all
clc
%% LOADING THE DATA
if ~exist('BZ_tensor')
    load('BZ.mat');
end


[m,n,k]=size(BZ_tensor); % x vs y vs time data


%% RE-SAMPLING THE DATA
% Reducing z size of BZ_tensor by sampling the frames
Sampling_dt=10;

% Number of sampled frames ks
ks=floor(k/Sampling_dt);

%define a new state variable
BZ_tensor_s=zeros(m,n,ks);

%initialize the new state variable
for jjj=1:ks
    BZ_tensor_s(:,:,jjj) = BZ_tensor(:,:,jjj*Sampling_dt);
end

BZ_s = zeros(m*n,ks);
for jj=1:ks
    temp = BZ_tensor_s(:,:,jj);
    BZ_s(:,jj) = temp(:);
end


%% DEFINE HANKEL MATRIX
t_embed=60;

H = Hmatrix(BZ_s,t_embed);
kH=ks-t_embed;

% SVD of Hankel matrix
[U,S,V]=svd(H,'econ');

plot(diag(S)/sum(diag(S)),'o')
sgtitle(strcat('Singular',{' '},'Values',{' '},'Of',{' '},'Hankel',{' '},'Matrix'))
%% DMD on the Hankel Matrix

H1 = H(:,2:end);
H2 = H(:,1:end-1);

rank = 30;
dt = 1;

% Step 1
[U2,Sigma2,V2] = svd(H2,'econ');
U = U2(:,1:rank);
V = V2(:,1:rank);
Sigma = Sigma2(1:rank,1:rank);

% Step 2
A_tilde = U'*H1*V/Sigma;

% Step 3
[W,D] = eig(A_tilde);

% Step 4
Phi = H1*V/Sigma*W;
mu = diag(D);
omega = log(mu)/dt;


plot(real(diag(Sigma))./sum(real(diag(Sigma))),'o')
%%
figure()
plot(real(diag(D)),imag(diag(D)),'o')
sgtitle(strcat('Eigenvalues',{' '},'Of',{' '},'A-Tilde'))
xlabel('Real'),ylabel('Imaginary')
%%
%%Project data in the low-rank subspace
H_tilde=U'*H;




%% DYNAMICS RECONSTRUCTION/PREDICTION BASED ON LOW RANK TRUNCATION
H_tilde_future=zeros(size(H_tilde));

H_tilde_future(:,1)=A_tilde*H_tilde(:,1);
for i=2:ks
    H_tilde_future(:,i)=A_tilde*H_tilde_future(:,i-1);
end



%% BACK PROJECTION TO THE ORIGINAL DATA SPACE
H_future=U*H_tilde_future;
BZ_s_future=H_future((1:m*n),:);

BZ_s_tensor_future=reshape(BZ_s_future,[m,n,ks]);

%% COMPARISON - first half of the time duration

for j=1:kH
A=BZ_tensor_s(:,:,j);
subplot(1,2,1)
 pcolor(A), shading interp, pause(0.01)
B=BZ_s_tensor_future(:,:,j);
subplot(1,2,2)
sgtitle(strcat('frame',{' '},num2str(j)))
pcolor(B), shading interp, pause(0.01)
end
%% COMPARISON - second half of the time duration

figure()
for j=kH:ks
A=BZ_tensor_s(:,:,j);
subplot(1,2,1)
 pcolor(A), shading interp, pause(0.01)
B=BZ_s_tensor_future(:,:,j);
subplot(1,2,2)
sgtitle(strcat('frame',{' '},num2str(j)))
pcolor(B), shading interp, pause(0.01)
end




