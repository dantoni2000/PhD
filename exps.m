clear all
close all

n = 100;
l = 50;
p = 0;
c12 = 0; c13 = 0; c23 = 0;
sigma = 1e-6; 

g = linspace(1,l,l)';



G1 = diag([ones(1,l), sigma*ones(1,n-l)]);
invLambda11 = diag(1./diag(G1(1:l,1:l)));
Lambda21 = G1(l+1:end,l+1:end);

% G1 = diag([ones(1,l), sigma, zeros(1,n-l-1)]);
% invLambda11 = diag(1./diag(G1(1:l,1:l)));
% Lambda21 = G1(l+1:end,l+1:end);


g2 = (-g'+g(l,1)+1)./g(l,1);
G2 = diag([max(g2,sigma), sigma*ones(1,n-l)]);
invLambda12 = diag(1./diag(G2(1:l,1:l)));
Lambda22 = G2(l+1:end,l+1:end);


g3 = (1./(g').^5);
G3 = diag([max(g3,sigma), sigma*ones(1,n-l)]);
invLambda13 = diag(1./diag(G3(1:l,1:l)));
Lambda23 = G3(l+1:end,l+1:end);


g4 = (1./(g').^.5);
G4 = diag([max(g4,sigma), sigma*ones(1,n-l)]);
invLambda14 = diag(1./diag(G4(1:l,1:l)));
Lambda24 = G4(l+1:end,l+1:end);

for tries = 1:100

    exp1A = 0; exp2A = 0; exp3A = 0; exp4A = 0;
    exp1B = 0; exp2B = 0; exp3B = 0; exp4B = 0;
    PSz1 = 0; PSz2 = 0; PSz3 = 0; PSz4 = 0;
    sPSz1 = 0; sPSz2 = 0; sPSz3 = 0; sPSz4 = 0;
    
    % Espressione vera
    for t = 1:100
        Omega = randn(n,l+p);
        Omega1 = Omega(1:l,1:end);
        Omega2 = Omega(l+1:end,1:end);
        Z1 = Lambda21.^0.5 * Omega2 * pinv(Omega1) * invLambda11.^0.5;
        S1 = Lambda21.^0.5 * Omega2 * pinv(Omega1);
        PSz1 = PSz1 + norm((S1'*S1)/(eye(l)+Z1'*Z1), 'fro');
        
        Z2 = Lambda22.^0.5 * Omega2 * pinv(Omega1) * invLambda12.^0.5;
        S2 = Lambda22.^0.5 * Omega2 * pinv(Omega1);
        PSz2 = PSz2 + norm((S2'*S2)/(eye(l)+Z2'*Z2),'fro');
        
        Z3 = Lambda23.^0.5 * Omega2 * pinv(Omega1) * invLambda13.^0.5;
        S3 = Lambda23.^0.5 * Omega2 * pinv(Omega1);
        PSz3 = PSz3 + norm((S3'*S3)/(eye(l)+Z3'*Z3),'fro');

        Z4 = Lambda24.^0.5 * Omega2 * pinv(Omega1) * invLambda14.^0.5;
        S4 = Lambda24.^0.5 * Omega2 * pinv(Omega1);
        PSz4 = PSz4 + norm((S4'*S4)/(eye(l)+Z4'*Z4),'fro');
     
        sZ1 = Lambda21.^0.5 * Omega2 * pinv(Omega1) * invLambda11.^0.5;
        sS1 = Lambda21.^0.5 * Omega2 * pinv(Omega1);
        sPSz1 = sPSz1 + norm((sS1'*sS1)/(1+norm(sZ1'*sZ1)),'fro');
        
        sZ2 = Lambda22.^0.5 * Omega2 * pinv(Omega1) * invLambda12.^0.5;
        sS2 = Lambda22.^0.5 * Omega2 * pinv(Omega1);
        sPSz2 = sPSz2 + norm((sS2'*sS2)/(1+norm(sZ2'*sZ2)),'fro');
        
        sZ3 = Lambda23.^0.5 * Omega2 * pinv(Omega1) * invLambda13.^0.5;
        sS3 = Lambda23.^0.5 * Omega2 * pinv(Omega1);
        sPSz3 = sPSz3 + norm((sS3'*sS3)/(1+norm(sZ3'*sZ3)),'fro');
        % [~,s,~] = svd(sS3'*sS3);
        % sPSz3 = sPSz3 + norm(eye(k)*s(end,end)/(1+s(end,end)),'fro');

        sZ4 = Lambda24.^0.5 * Omega2 * pinv(Omega1) * invLambda14.^0.5;
        sS4 = Lambda24.^0.5 * Omega2 * pinv(Omega1);
        sPSz4 = sPSz4 + norm((sS4'*sS4)/(1+norm(sZ4'*sZ4)),'fro');
    end


    PSz1 = 1/100 * PSz1; PSz2 = 1/100 * PSz2; PSz3 = 1/100 * PSz3; PSz4 = 1/100 * PSz4;
    sPSz1 = 1/100 * sPSz1; sPSz2 = 1/100 * sPSz2; sPSz3 = 1/100 * sPSz3; sPSz4 = 1/100 * sPSz4;


    ratio1 = PSz1 / sPSz1;
    ratio2 = PSz2 / sPSz2;
    ratio3 = PSz3 / sPSz3;
    ratio4 = PSz4 / sPSz4;
    r1(tries) = ratio1;
    r2(tries) = ratio2;
    r3(tries) = ratio3;
    r4(tries) = ratio4;

end

semilogy(r1)
xlabel('try')
ylabel('rapporto costanti empiriche C1/C2')
legend('C1/C2')

figure
semilogy(r2)
xlabel('try')
ylabel('rapporto costanti empiriche C1/C3')
legend('C1/C3')

figure
semilogy(r3)
xlabel('try')
ylabel('rapporto costanti empiriche C2/C3')
legend('C2/C3')

figure
semilogy(r4)
xlabel('try')
ylabel('rapporto costanti empiriche C2/C3')
legend('C2/C3')


figure
semilogy(r1)
hold on
semilogy(r2)
hold on
semilogy(r3)
hold on
semilogy(r4)
xlabel('try')
ylabel('Costante empirica C')
legend('C1', 'C2', 'C3', 'C4')

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
hold on
semilogy(r4)
xlabel('try')
ylabel('Costante empirica C')
legend('C1', 'C2', 'C3', 'C4')