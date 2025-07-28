clear all
close all

n = 100;
k = 20;
p = 0;
c12 = 0; c13 = 0; c23 = 0;
g = linspace(1,n,n)';

G1 = diag([ones(1,k), 1e-3*ones(1,n-k)]);
invLambda11 = diag(1./diag(G1(1:k,1:k)));
Lambda21 = G1(k+1:end,k+1:end);

% G2 = diag(1./sqrt(linspace(1,n,n)'));
G2 = diag((-g+g(n,1)+1)./g(n,1));
invLambda12 = diag(1./diag(G2(1:k,1:k)));
Lambda22 = G2(k+1:end,k+1:end);

A3 = zeros(n,n);
l = 1e-2; h = 1e+2; q = 1e+0; s = 1e-6;

for j=1:20
    x = sprand(n,1,0.01);
    A3 = A3 + (h/j^2*x) * x';
end
for j=21:40
    x = sprand(n,1,0.01);
    A3 = A3 + (q/j^2*x) * x';
end
for j=41:60
    x = sprand(n,1,0.01);
    A3 = A3 + (l/j^2*x) * x';
end
for j=61:80
    x = sprand(n,1,0.01);
    A3 = A3 + (s/j^2*x) * x';
end
[~,G3,~]=svd(A3);
invLambda13 = diag(1./diag(G3(1:k,1:k)));
Lambda23 = G3(k+1:end,k+1:end);

for tries = 1:100

    exp1A = 0; exp2A = 0; exp3A = 0; exp1B = 0; exp2B = 0; exp3B = 0; PSz1 = 0; PSz2 = 0; PSz3 = 0; sPSz1 = 0; sPSz2 = 0; sPSz3 = 0;
    
    % Espressione vera
    for t = 1:100
        Omega = randn(n,k+p);
        Omega1 = Omega(1:k,1:end);
        Omega2 = Omega(k+1:end,1:end);
        Z1 = Lambda21.^0.5 * Omega2 * pinv(Omega1) * invLambda11.^0.5;
        S1 = Lambda21.^0.5 * Omega2 * pinv(Omega1);
        PSz1 = PSz1 + norm((S1'*S1)/(eye(k)+Z1'*Z1), 'fro');
        Z2 = Lambda22.^0.5 * Omega2 * pinv(Omega1) * invLambda12.^0.5;
        S2 = Lambda22.^0.5 * Omega2 * pinv(Omega1);
        PSz2 = PSz2 + norm((S2'*S2)/(eye(k)+Z2'*Z2),'fro');
        Z3 = Lambda23.^0.5 * Omega2 * pinv(Omega1) * invLambda13.^0.5;
        S3 = Lambda23.^0.5 * Omega2 * pinv(Omega1);
        PSz3 = PSz3 + norm((S3'*S3)/(eye(k)+Z3'*Z3),'fro');

        sZ1 = Lambda21.^0.5 * Omega2 * pinv(Omega1) * invLambda11.^0.5;
        sS1 = Lambda21.^0.5 * Omega2 * pinv(Omega1);
        sPSz1 = sPSz1 + norm((sS1'*sS1)/(1+norm(sZ1'*sZ1)),'fro');
        sZ2 = Lambda22.^0.5 * Omega2 * pinv(Omega1) * invLambda12.^0.5;
        sS2 = Lambda22.^0.5 * Omega2 * pinv(Omega1);
        sPSz2 = sPSz2 + norm((sS2'*sS2)/(1+norm(sZ2'*sZ2)),'fro');
        sZ3 = Lambda23.^0.5 * Omega2 * pinv(Omega1) * invLambda13.^0.5;
        sS3 = Lambda23.^0.5 * Omega2 * pinv(Omega1);
        sPSz3 = sPSz3 + norm((sS3'*sS3)/(1+norm(sZ3'*sZ3)),'fro');
    end


    PSz1 = 1/100 * PSz1; PSz2 = 1/100 * PSz2; PSz3 = 1/100 * PSz3;
    sPSz1 = 1/100 * sPSz1; sPSz2 = 1/100 * sPSz2; sPSz3 = 1/100 * sPSz3;


    ratio1 = PSz1 / sPSz1;
    ratio2 = PSz2 / sPSz2;
    ratio3 = PSz3 / sPSz3;
    r1(tries) = ratio1;
    r2(tries) = ratio2;
    r3(tries) = ratio3;

    if ratio1/ratio2<100 && ratio2/ratio1<100
        c12 = c12+1;
    end

    if ratio1/ratio3<100 && ratio3/ratio1<100
        c13 = c13+1;
    end
    if ratio2/ratio3<100 && ratio3/ratio2<100
        c23 = c23+1;
    end
end

semilogy(r1)
xlabel('try')
ylabel('rapporto costanti empiriche C1/C2')
legend('C1/C2')

figure
semilogy(r3)
xlabel('try')
ylabel('rapporto costanti empiriche C1/C3')
legend('C1/C3')

figure
semilogy(r2)
xlabel('try')
ylabel('rapporto costanti empiriche C2/C3')
legend('C2/C3')

figure
semilogy(r1)
hold on
semilogy(r2)
hold on
semilogy(r3)
xlabel('try')
ylabel('Costante empirica C')
legend('C1', 'C2', 'C3')

figure
subplot('Position', [0.05 0.55 0.4 0.4])  % in alto a sinistra
semilogy(diag(G1))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

% Subplot 2: Nuclear norm error
subplot('Position', [0.55 0.55 0.4 0.4])  % in basso a destra
semilogy(diag(G2))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

% Subplot 3: Frobenius norm error
subplot('Position', [0.05 0.05 0.4 0.4])  % in alto a destra
semilogy(diag(G3))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

% Subplot 4: comparison
subplot('Position', [0.55 0.05 0.4 0.4])  % in alto a destra
semilogy(r1)
hold on
semilogy(r2)
hold on
semilogy(r3)
xlabel('try')
ylabel('Costante empirica C')
legend('C1', 'C2', 'C3')