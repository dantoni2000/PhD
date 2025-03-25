% 20 03 2025 Test Nystrom: 
% commenti stupidi: nei casi f operator monotone fcts prendo direttamente
% ANy, qui non posso perchè va a -infty. Perciò prendo ANy + I,
% direttamente facendo log det.

clear all

T = 10;
n = 1000;
Q = randn(n,n); [Q,~] = qr(Q);
g=linspace(1,n,n);
%G = diag(1./(g).^1);
G = diag(1./sqrt(g));
% G = diag(exp(-g));
A = Q*G*Q'; mi = 1e-10;
trA = sum(log(diag(G+mi*eye(n,n))),"all");
l = 20;

for tries=1:T
    for q=1:1:100
        
        [U,Lhat] = Nystrom(A,l+q);
        U = U(:,1:l); Lhat = Lhat(1:l,1:l); 
        eA(q) = trace(A-U*Lhat*U')/sum(diag(G(l+1:end,l+1:end))) - 1;
    
        % Ainv = Q*(Q'./diag(G+mi*eye(n,n))) + 1/mi*(eye(n,n)-Q*Q');
        % ANyinv = U*(U'./diag(Lhat+mi*eye(l,l))) + 1/mi*(eye(n,n)-U*U');
    
        trP(q) = abs(trA - log(prod(diag(Lhat+mi*eye(l,l)))))/(abs(trA - sum(log(diag(G(1:l,1:l)+mi*eye(l,l))),"all"))) - 1;
        
    
    
    end
    
    plot(eA,'g');
    hold on 
    plot(trP,'r');
    hold on
    legend('$\frac{||A - A_{2p-1}^{NYS}||_{*}}{||A - A_{2p-1}||_{*}}$','$\frac{||\log(A+ \mu I) - \log(A_{2p-1}^{NYS} + \mu I)||_{*}}{||\log(A + \mu I) -\log(A_{2p-1} + \mu I)||_{*}}$','interpreter','latex')
    % legend('$||A - A_{2p-1}^{NYS}||$','$||\lambda_{p}||$','interpreter','latex')
    % figure
    % loglog(eAinv)
    % legend('$||(A+\mu I)^{-1} - (A_{2p-1}^{NYS} +\mu I)^{-1}||$','interpreter','latex')
    % figure
    % loglog(AdeAinv)
    % legend('$||(A+\mu I)^{-1} - (A_{2p-1}^{AdNYS} +\mu I)^{-1}||$','interpreter','latex')
end