clear all
close all
warning off

T = 30;

% n = 500;
n = 1000;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n)';
mi = 1;
sigma = 1e-6;
nTr = zeros(21,1);
rk = 490;
G = diag([ones(1,rk), sigma*ones(1,n-rk)]);

% g = linspace(1,rk,rk)';
% g4 = (1./(g').^.5);
% G = diag([max(g4,sigma), sigma*ones(1,n-rk)]);

A = Q*G*Q'; 

for l=1:1:21
    ll = 479 + l;
    mvecs( l ) = ll;      
    for t=1:T
    % Nystrom grande su A
    [UBig,LhatBig] = PinvNystrom(A,ll);
    B = A-UBig*LhatBig*UBig';
    nTr(l,t) = trace(B);
    end
    

end
nTr = median(nTr,2);

for tMV = 1:1:21
    ttMV = 479 + tMV;
    BestRkTr(tMV) = sum(diag(G(ttMV+1:n,ttMV+1:n)));
    denom(tMV) = min( 1/ttMV + (3/G(ttMV,ttMV)) * BestRkTr(tMV), 1) ;
end

figure(1)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(2)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(2)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nTr'./BestRkTr,'-db')
hold on
semilogy(mvecs,1./denom,'-dg')
hold on
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for Nystrom in Trace Norm for $\sigma = 10^{-6}$', 'Interpreter','latex')
legend('$\frac{\|\Lambda-\Lambda_{k,p}\|_*}{\|\Lambda-\Lambda_{k+p}\|_*}$', '$\frac{k+p}{1+3(k+p)\|\sigma I \|_*}$', 'Fontsize', 25, 'interpreter', 'Latex')