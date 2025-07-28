clear all
close all

n = 100;
k = 20;
p = 0;
c12 = 0; c13 = 0; c23 = 0;
g = linspace(1,n,n)';

G1 = diag(1./(g).^5);
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

    exp1A = 0; exp2A = 0; exp3A = 0; exp1B = 0; exp2B = 0; exp3B = 0; PSz1 = 0; PSz2 = 0; PSz3 = 0;
    
    % primo termine
    for t = 1:100
        Omega = randn(n,k+p);
        Omega1 = Omega(1:k,1:end);
        Omega2 = Omega(k+1:end,1:end);
        Z1 = Lambda21.^0.5 * Omega2 * pinv(Omega1) * invLambda11.^0.5;
        exp1A = exp1A + 1/(1+norm(Z1'*Z1));
        Z2 = Lambda22.^0.5 * Omega2 * pinv(Omega1) * invLambda12.^0.5;
        exp2A = exp2A + 1/(1+norm(Z2'*Z2));
        Z3 = Lambda23.^0.5 * Omega2 * pinv(Omega1) * invLambda13.^0.5;
        exp3A = exp3A + 1/(1+norm(Z3'*Z3));
    end

    exp1A = 1/100 * exp1A; exp2A = 1/100 * exp2A; exp3A = 1/100 * exp3A;
    
    % secondo termine
    for t = 1:100
        Omega = randn(n,k+p);
        Omega1 = Omega(1:k,1:end);
        Omega2 = Omega(k+1:end,1:end);
        S1 = Lambda21.^0.5 * Omega2 * pinv(Omega1);
        exp1B = exp1B + norm(S1'*S1,'fro');
        S2 = Lambda22.^0.5 * Omega2 * pinv(Omega1);
        exp2B = exp2B + norm(S2'*S2,'fro');
        S3 = Lambda23.^0.5 * Omega2 * pinv(Omega1);
        exp3B = exp3B + norm(S3'*S3,'fro');
    end

    exp1B = 1/100 * exp1B; exp2B = 1/100 * exp2B; exp3B = 1/100 * exp3B;
    
    % Polya Szego
    for t = 1:100
        Omega = randn(n,k+p);
        Omega1 = Omega(1:k,1:end);
        Omega2 = Omega(k+1:end,1:end);
        Z1 = Lambda21.^0.5 * Omega2 * pinv(Omega1) * invLambda11.^0.5;
        S1 = Lambda21.^0.5 * Omega2 * pinv(Omega1);
        PSz1 = norm((S1'*S1)/(1+norm(Z1'*Z1)),'fro');
        Z2 = Lambda22.^0.5 * Omega2 * pinv(Omega1) * invLambda12.^0.5;
        S2 = Lambda22.^0.5 * Omega2 * pinv(Omega1);
        PSz2 = norm((S2'*S2)/(1+norm(Z2'*Z2)),'fro');
        Z3 = Lambda23.^0.5 * Omega2 * pinv(Omega1) * invLambda13.^0.5;
        S3 = Lambda23.^0.5 * Omega2 * pinv(Omega1);
        PSz3 = norm((S3'*S3)/(1+norm(Z3'*Z3)),'fro');
    end
    
    PSz1 = 1/100 * PSz1; PSz2 = 1/100 * PSz2; PSz3 = 1/100 * PSz3;

    ratio1 = PSz1 / (exp1A*exp1B);
    ratio2 = PSz2 / (exp2A*exp2B);
    ratio3 = PSz3 / (exp3A*exp3B);
    r12(tries) = ratio1/ratio2;
    r13(tries) = ratio1/ratio3;
    r23(tries) = ratio2/ratio3;
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

semilogy(r12)
xlabel('try')
ylabel('rapporto costanti empiriche C1/C2')
legend('C1/C2')

figure
semilogy(r13)
xlabel('try')
ylabel('rapporto costanti empiriche C1/C3')
legend('C1/C3')

figure
semilogy(r23)
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