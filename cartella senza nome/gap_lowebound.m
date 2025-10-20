clear all
close all
warning off

T = 30;

% n = 500;
n = 1500;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n)';
mi = 1;
sigma = 1e-5;
nTr = zeros(20,1);
rk = 290;
G = diag([ones(1,rk), sigma*ones(1,rk), sigma^2*ones(1,rk), sigma^4*ones(1,n-3*rk)]);

% g = linspace(1,rk,rk)';
% g4 = (1./(g').^.5);
% G = diag([max(g4,sigma), sigma*ones(1,n-rk)]);

A = Q*G*Q'; 

j = [289 290 291 292 300 579 580 581 582 590 869 870 871 872 880 1159 1160 1161 1162 1170]; 
for l=1:1:20
    ll = j(l);
    mvecs( l ) = ll;      
    for t=1:T
    % Nystrom grande su A
    [UBig,LhatBig] = PinvNystrom(A,ll);
    B = A-UBig*LhatBig*UBig';
    nTr(l,t) = trace(B);
    end
    

end
nTr = median(nTr,2);

for tMV = 1:1:20
    ttMV = j(tMV);
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
semilogy(diag(G), 'LineWidth', 5)
xlabel('$n$','fontsize',18,'interpreter','Latex')
ylabel('eigenvalues','fontsize',18)
title('Eigenvalues of the matrix','fontsize',18)
legend('$\lambda(A)$','fontsize',18,'interpreter','Latex')

figure(2)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs(1,1:20),nTr(1:20,1)'./BestRkTr(1,1:20),'-db', 'LineWidth', 5)
hold on
semilogy(mvecs(1,1:20),1./denom(1,1:20),'-dg', 'LineWidth', 5)
hold on
xlabel('MatVecs','fontsize',18)
ylabel('error','fontsize',18)
title('Comparison of the bounds for Nystrom in Trace Norm for $\sigma = 10^{-6}$','fontsize',18, 'Interpreter','latex')
legend('$\frac{\|\Lambda-\Lambda_{k,p}\|_*}{\|\Lambda-\Lambda_{k+p}\|_*}$', '$\frac{k+p}{1+3(k+p)\|\sigma I \|_*}$', 'fontsize',18, 'interpreter', 'Latex')