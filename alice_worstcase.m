clear all
close all
warning off

T = 50;

n = 500;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n)';
mi = 1;


nFro = zeros(49,1);
nTr = zeros(49,1);

for l=10:10:490

    mvecs( l/10 ) = l;
    G = diag([ones(1,l), 1e-3*ones(1,n-l)]);
    A = Q*G*Q';        
    for t=1:T
    % Nystrom grande su A
    [UBig,LhatBig] = PinvNystrom(A,l);
    B = A-UBig*LhatBig*UBig';
    nFro(l/10) = nFro(l/10) + norm(B,'fro');
    nTr(l/10) = nTr(l/10) + trace(B);
    end
    

end

nFro = 1/T * nFro;
nTr = 1/T * nTr;

for tMV = 10:10:490
    G = diag([ones(1,tMV), 1e-3*ones(1,n-tMV)]);
    BestRkTr(tMV/10) = sum(diag(G(tMV+1:n,tMV+1:n)));
    BestRkMENO1Tr = sum(diag(G(tMV:n,tMV:n)));
    BestRkFro(tMV/10) = sum(diag(G(tMV+1:n,tMV+1:n)).^2)^0.5;
    BestRkMENO1Fro(tMV/10) = sum(diag(G(tMV:n,tMV:n)).^2)^0.5;
    gap(tMV/10) = G(tMV,tMV)/G(tMV+1,tMV+1);
end

figure(1)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(1)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nFro'./BestRkFro,'-dr')
hold on
semilogy(mvecs,nFro'./BestRkMENO1Fro,'-dm')
hold on
semilogy(mvecs,gap,'-dk')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for Nystrom in Frobenius Norm')
legend('ratio error vs best rank k', 'ratio error vs best rank k-1', 'spectral gap')


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
semilogy(mvecs,nTr'./BestRkMENO1Tr,'-dc')
hold on
semilogy(mvecs,sqrt(gap),'-dk')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for Nystrom in Trace Norm')
legend('ratio error vs best rank k', 'ratio error vs best rank k-1', 'spectral gap')