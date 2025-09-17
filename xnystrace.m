function [t,err] = xnystrace(A, m)
N = size(A,1);

cnormc = @(M) M./ vecnorm(M,2,1);
diag_prod = @(A,B) sum(conj(A).*B,1).';

% [Q,B] = PinvNystrom(A,m);

Om = sqrt(N) * cnormc(randn(N,m));
Y = A*Om;
nu = eps*norm(Y,'fro')/sqrt(N);
Y = Y + nu*Om;
[Q,R] = qr(Y,0);
H = Om'*Y; C = chol((H+H')/2);
B = R/C;

[QQ,RR] = qr(Om,0);
WW = QQ'*Om;
SS = cnormc(inv(RR)');
scale = (N-m+1)./(N-vecnorm(WW).^2 + abs(diag_prod(SS,WW)'.*vecnorm(SS)).^2);

W = Q'*Om; S = (B/C').*(diag(inv(H))').^(-0.5);
dSW = diag_prod(S,W).';
ests = norm(B,'fro')^2 - vecnorm(S).^2 + abs(dSW).^2 .* scale - nu*N;
t = mean(ests);
err = std(ests)/sqrt(m);