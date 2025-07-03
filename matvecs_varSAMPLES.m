% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off

T = 30;

decadimento=7;

switch decadimento
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
        G = diag(exp(-0.1*g));
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
        G = sparse(diag(exp(-(1:n)/500)));
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
        [~,G]=eig(A); G = real(G);
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
        [~,G] = eig(A);
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

        case 10
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        % mi = 1;
        G = diag(1./(g).^.9);
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
end

trA = sum(log(diag(G+mi*eye(n,n))),"all");

ExperrtrBig=zeros(60,1);
Experrtr1=zeros(60,1);
ExperrtrNoPrec = zeros(60,1);

for t=1:T

    for l=10:10:600

        mvecs(l/10) = l;

        % Nystrom grande su A
        [UBig,LhatBig] = PinvNystrom(A,l);
        PBig = (mi)^0.5*UBig*(LhatBig+mi*eye(l))^-0.5*UBig' + (eye(n) - UBig*UBig');
        trBig(l/10,t) = sum(log(diag(LhatBig+mi*eye(l,l))),"all");
        ExperrtrBig(l/10) = ExperrtrBig(l/10) + abs(trBig(l/10,t) - trA);
        % trBig(l - 9 ) = log(prod((llBig+mi)*diag(LhatBig+mi*eye(l,l))));
        % trBig(l - 9 -1,t) = -2*trace(logm((llBig+mi)^-0.5*PBig));
    
        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-5;
        [Nits(l/10),tr1(l/10,t)] = Nystrom_HUTCH(A,mi,l1,1,50,1,1);
        Experrtr1(l/10) = Experrtr1(l/10) + abs(tr1(l/10,t) - trA);
        

        %no Nystrom, l/5 Hutch 5 Lanczos
        [its(l/10),tr_5HL] = EST(A,mi,l/5,50,1);
        trA_5HL(l/10,t) = tr_5HL;
        ExperrtrNoPrec(l/10) = ExperrtrNoPrec(l/10) + abs(trA_5HL(l/10,t) - trA);
        
    end

end

%mean errors

etrBig = 1/T * ExperrtrBig;
etr1 = 1/T * Experrtr1;
etrNoPrec = 1/T * ExperrtrNoPrec;

figure(3)
semilogy(mvecs,etrBig,'o-b')
hold on
semilogy(mvecs,etr1,'o-r')
hold on
semilogy(mvecs,etrNoPrec,'o-g')


% legend('Big Nystrom l', 'Nystrom l-25 5Hutch 5Lanczos', 'Nystrom l-5 1Hutch 5Lanczos', 'Nystrom l-10 2Hutch 5Lanczos', 'Nystrom l-30 3Hutch 10Lanczos', 'l - 9 samples 10 Lanczos', 'l/5samples 5 Lanczos') 
legend('Error Big Nystrom', 'Error Nystrom + Lanczos', 'Error No Nystrom') 
xlabel('matvecs')
ylabel('error')
title('error with respect to matvecs')

% % intervals
% Max_tr5 = max(tr5');
% Min_tr5 = min(tr5');
% 
% Max_trBig = max(trBig');
% Min_trBig = min(trBig');
% 
% Max_tr1 = max(tr1');
% Min_tr1 = min(tr1');
% 
% Max_tr2 = max(tr2');
% Min_tr2 = min(tr2');
% 
% Max_trA_L = max(trA_L');
% Min_trA_L = min(trA_L');
% 
% Max_trA_HL = max(trA_HL');
% Min_trA_HL = min(trA_HL');
% 
% figure(3)
% plot(mvecs,Max_tr5,'r')
% hold on
% RR = plot(mvecs,mtr5,'o-r')
% hold on
% plot(mvecs,Min_tr5,'r')
% hold on
% semilogy(mvecs,Max_trBig,'m')
% hold on
% MM = semilogy(mvecs,mtrBig,'o-m')
% hold on
% semilogy(mvecs,Min_trBig,'m')
% hold on
% semilogy(mvecs,Max_tr1,'g')
% hold on
% GG = semilogy(mvecs,mtr1,'o-g')
% hold on
% semilogy(mvecs,Min_tr1,'g')
% hold on
% semilogy(mvecs,Max_tr2,'k')
% hold on
% YY = semilogy(mvecs,mtr2,'o-k')
% hold on
% semilogy(mvecs,Min_tr2,'k')
% % hold on
% % plot(mvecs,Max_trA_L,'c')
% % hold on
% % CC = plot(mvecs,mtrA_L,'o-c')
% % hold on
% % plot(mvecs,Min_trA_L,'c')
% % hold on
% % plot(mvecs,Max_trA_HL,'b')
% % hold on
% % BB = plot(mvecs,mtrA_HL,'o-b')
% % hold on
% % plot(mvecs,Min_trA_HL,'b')
% % 
% legend([MM, GG, YY],'Big Nystrom', 'Nystrom 1Hutch 5Lanczos', 'Nystrom 1Hutch 3Lanczos')
% xlabel('matvecs')
% ylabel('traces')
% title('error wrt matvecs')
% % 
% % % figure(4)
% % % plot(mvecs,abs(Max_trBig - trA)./abs(trA),'r')
% % % hold on
% % % RR = plot(mvecs,abs(mtrBig - trA)./abs(trA),'o-r')
% % % hold on
% % % plot(mvecs,abs(Min_trBig - trA)./abs(trA),'r')
% % % hold on
% % % plot(mvecs,abs(Max_tr1 - trA)./abs(trA),'g')
% % % hold on
% % % GG = plot(mvecs,abs(mtr1 - trA)./abs(trA),'o-g')
% % % hold on
% % % plot(mvecs,abs(Min_tr1 - trA)./abs(trA),'g')
% % % hold on
% % % plot(mvecs,abs(Max_tr2 - trA)./abs(trA),'y')
% % % hold on
% % % YY = plot(mvecs, abs(mtr2 - trA)./abs(trA),'o-y')
% % % hold on
% % % plot(mvecs,abs(Min_tr2 - trA)./abs(trA),'y')
% % % hold on
% % % plot(mvecs,abs(Max_trA_L - trA)./abs(trA),'c')
% % % hold on
% % % CC = plot(mvecs,abs(mtrA_L - trA)./abs(trA),'o-c')
% % % hold on
% % % plot(mvecs,abs(Min_trA_L - trA)./abs(trA),'c')
% % % hold on
% % % plot(mvecs,abs(Max_trA_HL - trA)./abs(trA),'b')
% % % hold on
% % % BB = plot(mvecs,abs(mtrA_HL - trA)./abs(trA),'o-b')
% % % hold on
% % % plot(mvecs,abs(Min_trA_HL - trA)./abs(trA),'b')