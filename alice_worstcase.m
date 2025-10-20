clear all
close all
warning off

T = 100;

% n = 500;
n = 1000;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n)';
mi = 1;
sigma = 1e-6;
nTr = zeros(49,1);

for l=10:10:490

    mvecs( l/10 ) = l;
    G = diag([ones(1,l), sigma*ones(1,n-l)]);
    A = Q*G*Q';        
    for t=1:T
    % Nystrom grande su A
    [UBig,LhatBig] = PinvNystrom(A,l);
    B = A-UBig*LhatBig*UBig';
    nTr(l/10,t) = trace(B);
    end
    

end
nTr = median(nTr,2);

for tMV = 10:10:490
    G = diag([ones(1,tMV), sigma*ones(1,n-tMV)]);
    BestRkTr(tMV/10) = sum(diag(G(tMV+1:n,tMV+1:n)));
    denom(tMV/10) = 1/tMV + 3 * sigma *(n-tMV);
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
semilogy(diag(G), 'LineWidth', 5)
xlabel('$n$','fontsize',18,'interpreter','Latex')
ylabel('eigenvalues','fontsize',18)
title('Eigenvalues of the matrix','fontsize',18)
legend('$\lambda(A)$','fontsize',18,'interpreter','Latex')

figure(2)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nTr'./BestRkTr,'-db', 'LineWidth', 5)
hold on
semilogy(mvecs,1./denom,'-dg', 'LineWidth', 5)
hold on
xlabel('MatVecs','fontsize',18)
ylabel('error','fontsize',18)
title('Comparison of the bounds for Nystrom in Trace Norm for $\sigma = 10^{-6}$','fontsize',18, 'Interpreter','latex')
legend('$\frac{\|\Lambda-\Lambda_{k,p}\|_*}{\|\Lambda-\Lambda_{k+p}\|_*}$', '$\frac{k+p}{1+3(k+p)\|\sigma I \|_*}$', 'Fontsize', 18, 'interpreter', 'Latex')