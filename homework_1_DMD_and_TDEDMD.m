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

[m,L]=size(X)

%last snapshot
xk = X(:,L)

[X0, X1] = delay(X)

[Phi, Lambda, b] = DMD(X0,X1, 2)


%full DMD matrix A
A = X1*pinv(X0);

%prediction: next snapshot
% xkplus1 = A*xk
% xkplus2 = A*xkplus1
% xkplus3 = A*xkplus2

%long prediction

yearsplus=zeros(2*L,1);

for j= 1:2*L
    yearsplus(j)= 1845 + 2*(j-1);
end

Xplus=X
Xplus2=X

for j= L+1:2*L
    Xplus(:,j)= A*Xplus(:,j-1);
    Xplus2(:,j)= Phi*Lambda^(j-1)*b;
end
% figure()
% plot(Xplus')
% figure()
% plot(Xplus2')

populationDMD=real(Xplus');

figure()
plot(yearsplus, populationDMD(:,1),'-or', yearsplus, populationDMD(:,2),'-ob')
ylabel('Population [Thousands]') 
xlabel('Year') 
legend('Hare','Lynx')
title('Canadian lynx and snowshoe hare populations and prediction')


% eigenvalue decomposition of A for interpretability
[W,D] = eig(A);
%W eigenvectors
%D eigenvalues
figure()
plot(real(diag(D)),imag(diag(D)),'o')
sgtitle(strcat('Eigenvalues',{' '},'Of',{' '},'A'))
xlabel('Real'),ylabel('Imaginary')

%Regular DMD
%[u,s,v]=svd(x,'econ');
[u,s,v]=svd(X);
figure()
plot(diag(s)./sum(diag(s)),'ok')

%y=inv(u)*x;

figure()
plot(X(:,1),'-or')

%% TIME DELAY EMBEDDING DMD:

% Create Hankel matrix H:
H = Hmatrix(X,11);

% SVD of Hankel matrix H:
[Hu,Hs,Hv]=svd(H);

figure()
plot(diag(Hs),'or')

% Create H0, H1 1 step delayed matrices
[H0, H1] = delay(H);

% DMD of Henkel Matrix truncated to rank r
r=11;

[HPhi, HLambda, Hb] = DMD(H0, H1, r)

HPhiX=HPhi(1:2,:)


% TDE prediction
Xplus3=X

for j= L+1:2*L
    Xplus3(:,j)= HPhiX*HLambda^(j-1)*Hb;
end

populationTDDMD=real(Xplus3');

figure()
plot(yearsplus, populationTDDMD(:,1),'-or', yearsplus, populationTDDMD(:,2),'-ob')
ylabel('Population [Thousands]') 
xlabel('Year') 
legend('Hare','Lynx')
title('Canadian lynx and snowshoe hare populations and prediction')
