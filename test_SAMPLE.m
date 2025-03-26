%26 03 2025: test sul numero di sample e iterazioni di Lanczos per avere
%uno stimatore vicino epsilon a trace(log(A))

clear all

T=10;
n = 2000;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n);
% G = diag(1./(g).^1);
G = diag(1+1./sqrt(g));
% G = diag(1+exp(-g));
A = Q*G*Q'; mi = 1; %mi = 1e-10;
trA = sum(log(diag(G)),"all");

alpha = 1.1; ep = .1; delta = .1;
N = fix((4*alpha^2)*ep^-2*log(4/delta))+1; %manca il termine con la norma di log(A) % *(norm(log(diag(G)))^2+ep*norm(log(diag(G))))
m = fix(sqrt(G(1,1)/G(n,n)+1)/4*log(n^2/2*(alpha/(alpha-1))*sqrt(G(1,1)/G(n,n)+1)*log(2*(G(1,1)/G(n,n)))))+1;

for t=1:T
    tr(t) = 0;
    
    for j=1:N
        
        x=randn(n,1);
        tr(t) = tr(t) + Lanczos_log(A,x,m);
        
    end
    
    trN(t)= 1/N*tr(t);
    
    err_est(t)=abs(trN(t) - trA);
end

plot([1:T],err_est)