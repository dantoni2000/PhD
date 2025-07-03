clear all
close all
warning off
T = 20;

decadimento=1;

switch decadimento
    case 10
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = eye(n);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

    case 0
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = diag((-g+g(n,1)+1)./g(n,1));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

    case 1
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        % mi = 1;
        G = diag(1./(g).^5);
        A = Q*G*Q';
        % A = G;
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

    case 2
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = diag(1./sqrt(g));
        A = Q*G*Q';
        % A = G;
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

    case 3
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = diag(exp(-.1*g));
        A = Q*G*Q';
        % A = G;
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

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
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar
    
    case 5
        n = 1000;
        mi = 1;
        alpha = 1; nu = 5/2; % nu = 13/2;
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
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

    case 6
        n = 1000; mi = 1;
        G = 100*sparse(diag((1:n).^(-2)));
        % Q = gallery('orthog',n,1);
        [Q,~] = qr(randn(n,n));
        A = Q*G*Q';  
        % A = G;
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')

        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

    case 7
        n = 1000; mi = 1;
        G = sparse(diag(exp(-(1:n)/100)));
        % D = sparse(diag((1:matrix_size)).^(-3));
        % Q = gallery('orthog',n,1);
        [Q,~] = qr(randn(n,n));
        A = Q*G*Q'; 
        % A = G;
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

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
        figure(1)
        [~,G,~]=svd(A);
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

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
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        xlabel('$n$','interpreter','Latex')
        ylabel('eigenvalues')
        title('Eigenvalues of kernel matrix')
        legend('$\lambda(A)$','$\lambda(A + \mu I)$','interpreter','Latex')
                
        figure(2)
        imagesc(log(abs(A)), [-16,3])
        colorbar

end

nFro = zeros(19,1);
nTr = zeros(19,1);

for l=10:10:190

    mvecs( l/10 ) = l;
            
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

for tMV = 10:10:190

    for s = 4:2:tMV-4
        k = s;
        p = tMV - s;
        
        boundTr(:,s/2-1) = (1 + k/(p-1)).* norm(G(k+1:n,k+1:n).^0.5 ,'fro')^2;
        boundFro(:,s/2-1) = (1 + sqrt(k/(p-1))).^2.* norm(G(k+1:n,k+1:n),'fro') ;
        % 1/(1 + 1/sqrt(G(tMV,tMV)) *sqrt(G(tMV+1,tMV+1))*sqrt(tMV/(p-1)) + 1/sqrt(G(tMV,tMV)) *norm(diag(G(tMV+1:n,tMV+1:n)).^0.5,'fro').^2*exp(1)*sqrt(tMV+p)/(p-1))
        LowerboundTr(:,s/2-1) = norm(diag(G(tMV+p+1:n,tMV+p+1:n)).^0.5,'fro').^2 ;
        LowerboundFro(:,s/2-1) = norm(diag(G(tMV+p+1:n,tMV+p+1:n)),'fro') ;
    end
    
    BestTr(tMV/10) = min(boundTr');
    BestFro(tMV/10) = min(boundFro');
    BestRkTr(tMV/10) = norm(diag(G(tMV+1:n,tMV+1:n)).^0.5,'fro').^2;
    BestRkFro(tMV/10) = norm(diag(G(tMV+1:n,tMV+1:n)),'fro');
    BestLowerTr(tMV/10) = max(LowerboundTr');
    BestLowerFro(tMV/10) = max(LowerboundFro');
    
end

figure
semilogy(mvecs,nFro,'r')
hold on
semilogy(mvecs,nTr,'b')
hold on
% semilogy(mvecs,BestFro','-or')
% hold on
% semilogy(mvecs,BestTr','-ob')
% hold on
semilogy(mvecs,BestRkFro','-*r')
hold on
semilogy(mvecs,BestRkTr','-*b')
hold on
semilogy(mvecs,BestLowerFro','-dr')
hold on
semilogy(mvecs,BestLowerTr','-db')

xlabel('MatVecs')
ylabel('error')
legend('error Nystrom+Lanczos', 'error Big Nystrom', 'conjecture Nystrom+Lanczos', 'bound Big Nystrom')

% figure(100)
% plot(diag(A))

% tipo qualcosa che dipende dal gap ex 1/(1+gap*l) (1+l) ||A-A_k+p+l||
for jj = 1:20
    ratio(jj) = nFro(end)./norm(diag(G(tMV+1+jj:n,tMV+1+jj:n)),'fro');
end
figure(4)
plot(ratio)