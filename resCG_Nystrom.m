% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off

T = 1;

decadimento=0;

switch decadimento
    case 0
        n = 700;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10^5*diag((-g+g(n,1)+1)./g(n,1)).^15;
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

    case 1
        n = 700;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        % mi = 1;
        G = 10^5*diag(1./(g).^2);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
    case 2
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = diag(1./sqrt(g));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
    case 3
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = diag(exp(-g));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
    case 4
        n = 1000;
        mi = 1;
        kernel = @(x,y) ((norm(x)/n)^10 + (norm(y)/n)^10)^(0.1);
        data_matrix = orth(randn(n,n))*diag(1:n);
        A = build_kernel_matrix(data_matrix,kernel);
        [Q,G] = svd(A);
        G=diag(10^5*diag(G));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
    
    case 5
        n = 1000;
        mi = 1;
        alpha = 1; nu = 5/2;
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
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

    case 6
        n = 1000; mi = 1;
        G = 100*sparse(diag((1:n).^(-2)));
        Q = gallery('orthog',n,1);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

    case 7
        n = 1000; mi = 1;
        G = sparse(diag(exp(-(1:n)/100)));
        %D = sparse(diag((1:matrix_size)).^(-3));
        Q = gallery('orthog',n,1);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

    case 8
        n = 1000;
        A = zeros(n,n); mi = 1;
        l=1; h=10000;
        for j=1:60
            x = sprand(n,1,0.01);
            A = A + (h/j^2*x) * x';
        end
        for j=61:300
            x = sprand(n,1,0.01);
            A = A + (l/j^2*x) * x';
        end
        figure(1)
        [~,G]=eig(A); G = real(G);
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

    case 9
        n = 1000;
        A = zeros(n,n); mi = 1;
        l = 1e-0; h = 1e+4; k = 1e+5; p = 1e+1;

        for j=1:60
            x = sprand(n,1,0.01);
            A = A + (h/j^2*x) * x';
        end
        for j=61:100
            x = sprand(n,1,0.01);
            A = A + (k/j^2*x) * x';
        end
        for j=101:200
            x = sprand(n,1,0.01);
            A = A + (l/j^2*x) * x';
        end
        for j=201:300
            x = sprand(n,1,0.01);
            A = A + (p/j^2*x) * x';
        end
        figure(1)
        [~,G]=eig(A); G = real(G);
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

        
end
b = A*ones(n,1);

for t=1:T

    for l=11:5:701

        mvecs((l-6)/5 ) = l;

        % cg con l its
        [~,~,rcg((l-6)/5)]=pcg(A + mi*eye(n),b,1e-10,l);
    
        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-5;
        [U,Lhat] = Nystrom(A,l1);
        Ps = U*(Lhat+mi*eye(l1))^0.5*U' + (mi+0*Lhat(end,end))^0.5*(eye(n) - U*U');
        [~,~,rpcg1((l-6)/5)] = pcg(A + mi*eye(n),b,1e-10,5,Ps,Ps);
        
        %Nystrom con 2 Hutch 5 Lanczos
        l2 = l-5;
        [U2,Lhat2] = Nystrom(A,l2);
        P2s = U2*(Lhat2+mi*eye(l2))*U2' + (mi+0*Lhat2(end,end))*(eye(n) - U2*U2');
        [~,~,rpcg2((l-6)/5)] = pcg(A + mi*eye(n),b,1e-10,5,P2s);

    end

end

figure(2)
semilogy(mvecs,rcg,'o-m')
hold on
% semilogy(mvecs,abs(mtr5-trA)./abs(trA),'r')
% hold on
semilogy(mvecs,rpcg1,'o-g')
hold on
semilogy(mvecs,rpcg2,'o-y')
% hold on
% semilogy(mvecs,abs(mtr6-trA)./abs(trA),'k')
% hold on
% semilogy(mvecs,abs(mtrA_L-trA)./abs(trA),'c')
% hold on
% semilogy(mvecs,abs(mtrA_HL-trA)./abs(trA),'b')

% legend('Big Nystrom l', 'Nystrom l-25 5Hutch 5Lanczos', 'Nystrom l-5 1Hutch 5Lanczos', 'Nystrom l-10 2Hutch 5Lanczos', 'Nystrom l-30 3Hutch 10Lanczos', '(l-6)/5 samples 10 Lanczos', '(l-6)/5samples 5 Lanczos') 
legend('lCG', 'Nystrom l-5 5CG', 'Nystrom l-10 10CG') 
xlabel('matvecs')
ylabel('error')
title('relative error wrt matvecs')