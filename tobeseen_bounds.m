clear all
close all
warning off

T = 30;

decadimento=3;

switch decadimento

    case 0
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10000 * diag((-g+g(n,1)+1)./g(n,1));
        A = Q*G*Q';
                
        % figure(2)
        % imagesc(log(abs(A)), [-16,3])
        % colorbar

    case 1
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        % mi = 1;
        G = diag(1./(g).^5);
        A = Q*G*Q';
        % A = G;

    case 2
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10000 * diag(1./sqrt(g));
        A = Q*G*Q';
        % A = G;

    case 3
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10000 * diag(exp(-.1*g));
        A = Q*G*Q';
        % A = G;

    case 4
        n = 1000;
        mi = 1;
        kernel = @(x,y) ((norm(x)/n)^10 + (norm(y)/n)^10)^(0.1);
        data_matrix = orth(randn(n,n))*diag(1:n);
        A = build_kernel_matrix(data_matrix,kernel);
        [Q,G] = svd(A);
        G=diag(diag(G));
        A = Q*G*Q'; 
        % A = G;
        
    case 5
        n = 1000;
        mi = 1;
        alpha = 1; nu = 1/2; % nu = 13/2;
        kernel = @(x,y) sqrt(pi)*((alpha*norm(x-y))^(nu)*besselk(nu,alpha*norm(x-y)))/(2^(nu-1)*alpha^(2*nu)*gamma(nu+0.5));
        data_matrix = randn(1,n);
        
        for row = 1:n

            for column = 1:n

                if row == column
                
                    A(row,column) = sqrt(pi)*gamma(nu)/(gamma(nu + 1/2)*alpha^(2*nu));

                else

                    A(row,column) = kernel(data_matrix(:,row),data_matrix(:,column));

                end

            end
        
        end

        [Q,G] = svd(A);
        G=diag(diag(G));
        A = Q*G*Q'; 
        % A = G;
        
    case 6
        n = 1000; mi = 1;
        G = 100*sparse(diag((1:n).^(-2)));
        % Q = gallery('orthog',n,1);
        [Q,~] = qr(randn(n,n));
        A = Q*G*Q';  
        % A = G;
        
    case 7
        n = 1000; mi = 1;
        G = sparse(diag(exp(-(1:n)/100)));
        % D = sparse(diag((1:matrix_size)).^(-3));
        % Q = gallery('orthog',n,1);
        [Q,~] = qr(randn(n,n));
        A = Q*G*Q'; 
        % A = G;
        
    case 8
        n = 1000;
        A = zeros(n,n); mi = 1;
        l = 1e-2; h = 1e+2; k = 1e+0; p = 1e-6;

        for j=1:20
            x = sprand(n,1,0.01);
            A = A + (h/j^2*x) * x';
        end
        for j=21:40
            x = sprand(n,1,0.01);
            A = A + (k/j^2*x) * x';
        end
        for j=41:60
            x = sprand(n,1,0.01);
            A = A + (l/j^2*x) * x';
        end
        for j=61:300
            x = sprand(n,1,0.01);
            A = A + (p/j^2*x) * x';
        end
        [~,G,~]=svd(A);
        
        case 9
        n = 1000;
        x = randn(1,n);
        mi = 1;
        for i = 1:n
            for j = 1:n
                A(i,j) = exp(-(x(i)-x(j))^2 / (2)^2);
            end
        end
        [~,G,~] = svd(A);
        
        case 10
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = eye(n);
        A = Q*G*Q';
        
end

trA = sum(log(diag(G+mi*eye(n,n))),"all");
nTr = zeros(9,1);
nFro = zeros(9,1);
nNoPrec = zeros(9,1);

for l=15:10:95
        
    mvecs((l-5)/10 ) = l;
        
    for t=1:T
        % Nystrom grande su A
        [UBig,LhatBig] = Nystrom(A,l);
        nTr((l-5)/10) = nTr((l-5)/10) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA);
        
        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-5;
        [~,trr] = Nystrom_HUTCH(A,mi,l1,1,5,1,1);
        nFro((l-5)/10) = nFro((l-5)/10) + abs(trr - trA);
        
        %Nystrom con n Hutch
        [~,trrNoPrec] = EST(A,mi,l/5,5,1);
        nNoPrec((l-5)/10) = nNoPrec((l-5)/10) + abs(trrNoPrec - trA);

    end

end

nFro = 1/T * nFro;
nTr = 1/T * nTr;
nNoPrec = 1/T * nNoPrec;

for tMV = 10:10:90

    for s = 4:2:tMV-4
        k = s;
        p = tMV - s;
        
        boundFro(:,s/2-1) = 2*norm(logm(G(k+1:n,k+1:n)+eye(n-k)),'fro') + 2*sqrt(k/(p-1).* log(G(1,1)+1) *norm(logm(G(k+1:n,k+1:n)+eye(n-k)).^0.5 ,'fro')^2);
        bound2Fro(:,s/2-1) = sqrt(2)*( norm(G(k+1:n,k+1:n),'fro') + k/(p-1) *norm(G(k+1:n,k+1:n).^0.5 ,'fro')^2);
        % boundFro(:,s/2-1) = sqrt(2) *( norm(logm(G(k+1:n,k+1:n)+eye(n-k)),'fro') + k/(p-1) *norm(logm(G(k+1:n,k+1:n)+eye(n-k)).^0.5 ,'fro')^2);
        boundTr(:,s/2-1) = (1+ (k+5)/(p-1) )*norm(logm(G(k+5+1:n,k+5+1:n)+eye(n-k-5)).^0.5,'fro').^2;
        bound2Tr(:,s/2-1) = (1+ (k+5)/(p-1) )*norm(G(k+5+1:n,k+5+1:n).^0.5,'fro').^2;
    end
    
    BestTr(tMV/10) = min(boundTr');
    Best2Tr(tMV/10) = min(bound2Tr');
    BestFro(tMV/10) = min(boundFro');
    Best2Fro(tMV/10) = min(bound2Fro');
    
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
semilogy(mvecs,nFro,'r')
hold on
semilogy(mvecs,BestFro','-om')
hold on
semilogy(mvecs,Best2Fro','-*k')
% hold on
% semilogy(mvecs,BestRkTr','-*b')
% hold on
% semilogy(mvecs,BestLowerFro','-dr')
% hold on
% semilogy(mvecs,BestLowerTr','-db')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for log Nystrom in Frobenius Norm')
legend('error Nystrom', 'bound log Nystrom', 'bound Nystrom')


figure(2)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(2)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nTr,'b')
hold on
semilogy(mvecs,BestTr','-oc')
hold on
semilogy(mvecs,Best2Tr','-*k')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for log Nystrom in Nuclear norm')
legend('error Nystrom', 'square log Nystrom', 'bound Nystrom')


figure(3)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(3)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nFro,'r')
hold on
semilogy(mvecs,BestFro','-om')
hold on
semilogy(mvecs,Best2Fro','-*r')
hold on
semilogy(mvecs,nTr,'b')
hold on
semilogy(mvecs,BestTr','-oc')
hold on
semilogy(mvecs,Best2Tr','-*b')
% hold on
% semilogy(mvecs,BestRkTr','-*b')
% hold on
% semilogy(mvecs,BestLowerFro','-dr')
% hold on
% semilogy(mvecs,BestLowerTr','-db')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for logtrace')
legend('error Nystrom+Hutch', 'bound log Nystrom+Hutch', 'bound Nystrom+Hutch','error Nystrom', 'square log Nystrom', 'bound Nystrom')