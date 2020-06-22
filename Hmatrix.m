function [H] = Hmatrix(X,r)

[n,m] = size(X);

if m < 2
    error('The number of columns of X must be larger than 1')
elseif r >= m
    error('r must be smaller than the number of columns of X')
elseif r < 1
    error('r must be larger than 0')
end

H = zeros(n*(r+1),m-r);

for kk = 0:r
    ii = n*kk+1;
    jj = kk+1;
    H(ii:(ii+n-1),:) = X(:,jj:(jj+m-r-1));
end



