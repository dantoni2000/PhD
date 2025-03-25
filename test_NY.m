% 20 03 2025 Test Nystrom

n = 1000;
Q = randn(n,n); [Q,~] = qr(Q);
g=linspace(1,n,n);
% G = diag(1./(g).^2);
G = diag(exp(-g));
A = Q*G*Q';

for l=5:5:50
    
    [U,Lhat] = Nystrom(A,l);
    
    eA(l/5) = norm(A-U*Lhat*U');
    p=fix((l+1)/2);
    rA(l/5) = G(p,p);
    mi=1;
    
    Ainv = Q*(Q'./diag(G+mi*eye(n,n))) + 1/mi*(eye(n,n)-Q*Q');
    ANyinv = U*(U'./diag(Lhat+mi*eye(l,l))) + 1/mi*(eye(n,n)-U*U');
    
    eAinv(l/5) = norm(Ainv-ANyinv);
    
    l_0 = l;
    l_max = n/2;
    q = 5;
    tol_e = 1e-3;
    tol_r = 1e-3;
    [AdU,AdLhat] = Adaptive_Nystrom(A,l_0,l_max,q,tol_e,tol_r,mi);
    AdANyinv = AdU*(AdU'./diag(AdLhat+mi*eye(size(AdLhat)))) + 1/mi*(eye(n,n)-AdU*AdU');
    AdeAinv(l/5) = norm(Ainv-AdANyinv);

end

figure
plot(eA);
hold on 
plot(rA);
legend('$||A - A_{2p-1}^{NYS}||$','$||\lambda_{p}||$','interpreter','latex')
figure
loglog(eAinv)
legend('$||(A+\mu I)^{-1} - (A_{2p-1}^{NYS} +\mu I)^{-1}||$','interpreter','latex')
figure
loglog(AdeAinv)
legend('$||(A+\mu I)^{-1} - (A_{2p-1}^{AdNYS} +\mu I)^{-1}||$','interpreter','latex')
