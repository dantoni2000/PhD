function [U,Lhat] = Nystrom(A,l)

[n,~] = size(A);
Om = randn(n,l);
[Om,~] = qr(Om,0);
Y = A*Om;
nu = sqrt(n)*eps(norm(Y,2)); 
Y_nu = Y + nu*Om;
C = chol(Om'*Y_nu);
B = Y_nu/C;
[U,S,~] = svd(B,0);
% Lhat = max(0,S^2-nu*eye(m,m));
Lhat = max(0,S^2 - nu*eye(size(S)));