clear all
close all
warning off

T = 20;

% n = 500;
n = 500;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n)';
mi = 1;
sigma = 1e-6;

nTr1 = zeros(49,1);
nTr2 = zeros(49,1);
nTr3 = zeros(49,1);

for l=10:10:490

    mvecs( l/10 ) = l;
    G1 = diag([ones(1,l), sigma*ones(1,n-l)]);
    G2 = diag([ones(1,l-1), sigma*ones(1,n-l+1)]);
    % G3 = diag([ones(1,l+1), sigma*ones(1,n-l-1)]);
    % G3 = diag([ones(1,l-7), sigma*ones(1,n-l+7)]);
    G3 = diag([ones(1,l/2), sigma*ones(1,n-l/2)]);
    A1 = Q*G1*Q'; 
    A2 = Q*G2*Q';
    A3 = Q*G3*Q';

    for t=1:T
    % Nystrom grande su A
    [UBig,LhatBig] = PinvNystrom(A1,l);
    B1 = A1-UBig*LhatBig*UBig';
    nTr1(l/10) = nTr1(l/10) + trace(B1);

    [UBig,LhatBig] = PinvNystrom(A2,l);
    B2 = A2-UBig*LhatBig*UBig';
    nTr2(l/10) = nTr2(l/10) + trace(B2);

    [UBig,LhatBig] = PinvNystrom(A3,l);
    B3 = A3-UBig*LhatBig*UBig';
    nTr3(l/10) = nTr3(l/10) + trace(B3);
    end
    

end

nTr1 = 1/T * nTr1;
nTr2 = 1/T * nTr2;
nTr3 = 1/T * nTr3;

for tMV = 10:10:490
    G1 = diag([ones(1,tMV), sigma*ones(1,n-tMV)]);
    G2 = diag([ones(1,tMV-1), sigma*ones(1,n-tMV+1)]);
    % G3 = diag([ones(1,tMV+1), sigma*ones(1,n-tMV-1)]);
    % G3 = diag([ones(1,tMV-7), sigma*ones(1,n-tMV+7)]);
    G3 = diag([ones(1,tMV/2), sigma*ones(1,n-tMV/2)]);
    BestRkTr1(tMV/10) = sum(diag(G1(tMV+1:n,tMV+1:n)));
    BestRkTr2(tMV/10) = sum(diag(G2(tMV+1:n,tMV+1:n)));
    BestRkTr3(tMV/10) = sum(diag(G3(tMV+1:n,tMV+1:n)));
end

figure(1)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G1))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(1)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nTr1'./BestRkTr1,'-dr')
hold on
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for Nystrom in Frobenius Norm')
legend('ratio error vs best rank k')


figure(2)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G2))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(2)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nTr2'./BestRkTr2,'-db')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for Nystrom in Trace Norm')
legend('ratio error vs best rank k', 'ratio error vs best rank k-1', 'bound with $\|\Omega_1^\dagger\|_F = \sqrt{5k}$', 'bound with $\|\Omega_1^\dagger\|_F = k$', 'spectral gap', 'interpreter', 'Latex')


figure(3)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G3))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(3)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nTr3'./BestRkTr3,'-db')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for Nystrom in Trace Norm')
legend('ratio error vs best rank k', 'ratio error vs best rank k-1', 'bound with $\|\Omega_1^\dagger\|_F = \sqrt{5k}$', 'bound with $\|\Omega_1^\dagger\|_F = k$', 'spectral gap', 'interpreter', 'Latex')