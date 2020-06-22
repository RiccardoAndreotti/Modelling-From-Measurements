clearvars -except BZ_tensor BZ
close all
clc

if ~exist('BZ_tensor')
    load('BZ.mat');
end

[m,n,k]=size(BZ_tensor); % x vs y vs time data

BZ = zeros(m*n,k);
for jj=1:k
    temp = BZ_tensor(:,:,jj);
    BZ(:,jj) = temp(:);
end

[U,S,V]=svd(BZ,'econ');
%%
plot(diag(S)/sum(diag(S)),'o')
sgtitle(strcat('Singular',{' '},'Values',{' '},'Of',{' '},'BZ',{' '},'States'))
%% Standard DMD on the first half of the BZ states

frames2train=600; %max=k
X1 = BZ(:,2:frames2train);
X2 = BZ(:,1:frames2train-1);

rank = 30;
dt = 1;

% Step 1
[U2,Sigma2,V2] = svd(X2,'econ');
U = U2(:,1:rank);
V = V2(:,1:rank);
Sigma = Sigma2(1:rank,1:rank);

% Step 2
A_tilde = U'*X1*V/Sigma;

% Step 3
[W,D] = eig(A_tilde);

% Step 4
Phi = X1*V/Sigma*W;
mu = diag(D);
omega = log(mu)/dt;
y0 = Phi\BZ(:,1);

plot(real(diag(Sigma))./sum(real(diag(Sigma))),'o')
sgtitle(strcat('Rank-truncated',{' '},'Singular',{' '},'Values',{' '},'Of',{' '},'BZ',{' '},'States'))

% projecting data to a low rank subspace  
BZ_tilde=U'*BZ;
%% COMPARING RECONSTRUCTED DYNAMICS WITH ORIGINAL  
% forcasting the future dynamics from final state of the data
BZ_tilde_reconstruct=zeros(size(BZ_tilde));

BZ_tilde_reconstruct(:,1)=A_tilde*BZ_tilde(:,1);
for i=2:k
    BZ_tilde_reconstruct(:,i)=A_tilde*BZ_tilde_reconstruct(:,i-1);
end

% re-projecting forecasts to a full rank space
BZ_reconstruct=U*BZ_tilde_reconstruct;

% reshaping from columns to video frames 
BZ_tensor_reconstruct=reshape(BZ_reconstruct,[m,n,k]);

%% ANIMATION OF THE ORIGINAL AND RECONSTRUCTED DATA up to mid-time and forecast till the end

for j=1:3:k
A=BZ_tensor(:,:,j);
subplot(1,2,1)
pcolor(A), shading interp, pause(0.005)
B=BZ_tensor_reconstruct(:,:,j);
subplot(1,2,2)
sgtitle(num2str(j))
pcolor(B), shading interp, pause(0.005)
end
%% PLOT THE ORIGINAL AND RECONSTRUCTED DATA up to mid-time and forecast till the end
figure()
subplot(1,4,1)
pcolor(BZ_tensor(:,:,600)), shading interp
subplot(1,4,2)
pcolor(BZ_tensor(:,:,610)), shading interp
subplot(1,4,3)
pcolor(BZ_tensor(:,:,630)), shading interp
subplot(1,4,4)
pcolor(BZ_tensor(:,:,700)), shading interp

figure()
subplot(1,4,1)
pcolor(BZ_tensor_reconstruct(:,:,600)), shading interp
subplot(1,4,2)
pcolor(BZ_tensor_reconstruct(:,:,610)), shading interp
subplot(1,4,3)
pcolor(BZ_tensor_reconstruct(:,:,630)), shading interp
subplot(1,4,4)
pcolor(BZ_tensor_reconstruct(:,:,700)), shading interp

%sgtitle(strcat('DATA',{' '},'vs',{' '},'DMD-FORECAST'))
