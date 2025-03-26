% 26 03 2025 Test Nystrom: come fare per evitare di fare Lanczos ma avere
% un risultato più accurato per shift piccoli.

clear all

T = 10;
n = 100;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n);
G = diag(1./(g).^2);
%G = diag(1./sqrt(g));
%G = diag(exp(-g));
A = Q*G*Q'; mi = 1e+2;
trA = sum(log(diag(G+mi*eye(n,n))),"all");
l = 40;
[U,Lhat] = Nystrom(A,l); ll=Lhat(l,l);
ANyinv = U*(U'./diag(Lhat+mi*eye(l,l))) + 1/mi*(eye(n,n)-U*U');

Est = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l)))) + n*log(mi);