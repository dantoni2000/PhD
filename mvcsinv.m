% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off

T = 50;

decadimento=0;

switch decadimento
    case 0
        n = 200;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-4;
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

    case 1
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-3;
        G = diag(1./(g).^2);
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
        mi = 1e-2;
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
        n = 200;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-3;
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
        mi = 1e-4;
        kernel = @(x,y) ((norm(x)/n)^10 + (norm(y)/n)^10)^(0.1);
        data_matrix = orth(randn(n,n))*diag(1:n);
        A = build_kernel_matrix(data_matrix,kernel);
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
    case 5
        n = 1000;
        mi = 1e-5;
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
end

trA = sum(log(diag(G+mi*eye(n,n))),"all");

for t=1:T

    for l=20:10:100

        mvecs(l/10-1) = l;

    
        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-5;
        [U1,Lhat1] = Nystrom(A,l1);  ll1 = Lhat1(l1,l1);
        P1 = (ll1+mi)^0.5*U1*(Lhat1+mi*eye(l1))^-0.5*U1' + (eye(n) - U1*U1');
        % trP1 = log((ll1+mi)^(-l1)*prod(diag(Lhat1+mi*eye(l1,l1))));
        trP1 = -2*trace(logm(P1));

        [~,Ptr1] = Preconditioned_HUTCH(A,mi,P1,1,5,1,1);
        tr1(l/10-1,t) = trP1 + Ptr1;
        
        [iU1,iLhat1] = Nystrom(inv(A+1e-15*eye(n)),l1); iLhat1 = inv(iLhat1); ill1 = iLhat1(l1,l1);
        iP1 = (ill1+mi)^0.5*iU1*(iLhat1+mi*eye(l1))^-0.5*iU1' + (eye(n) - iU1*iU1');
        % trP1 = log((ll1+mi)^(-l1)*prod(diag(Lhat1+mi*eye(l1,l1))));
        triP1 = -2*trace(logm(iP1));

        [~,iPtr1] = Preconditioned_HUTCH(A,mi,iP1,1,5,1,1);
        itr1(l/10-1,t) = triP1 + iPtr1;

        %Nystrom con 2 Hutch 5 Lanczos
        l2 = l-10;
        [U2,Lhat2] = Nystrom(A,l2);  ll2 = Lhat2(l2,l2);
        P2 = (ll2+mi)^0.5*U2*(Lhat2+mi*eye(l2))^-0.5*U2' + (eye(n) - U2*U2');
        % trP2 = log((ll2+mi)^(-l2)*prod(diag(Lhat2+mi*eye(l2,l2))));
        trP2 = -2*trace(logm(P2));

        [~,Ptr2] = Preconditioned_HUTCH(A,mi,P2,2,5,1,1);
        tr2(l/10-1,t) = trP2 + Ptr2;

        [iU2,iLhat2] = Nystrom(inv(A+1e-15*eye(n)),l2); iLhat2 = inv(iLhat2); ill2 = iLhat2(l2,l2);
        iP2 = (ill2+mi)^0.5*iU2*(iLhat2+mi*eye(l2))^-0.5*iU2' + (eye(n) - iU2*iU2');
        % trP1 = log((ll1+mi)^(-l1)*prod(diag(Lhat1+mi*eye(l1,l1))));
        triP2 = -2*trace(logm(iP2));

        [~,iPtr2] = Preconditioned_HUTCH(A,mi,iP2,2,5,1,1);
        itr2(l/10-1,t) = triP2 + iPtr2;
        
        %Nystrom con 5 Hutch 2 Lanczos
        [~,Ptr5] = Preconditioned_HUTCH(A,mi,P2,5,2,1,1);
        tr5(l/10-1,t) = trP2 + Ptr5;

        [~,iPtr5] = Preconditioned_HUTCH(A,mi,iP2,5,2,1,1);
        itr5(l/10-1,t) = triP2 + iPtr5;

        %no Nystrom, 1 Hutch l Lanczos
        [~,tr_L] = EST(A,mi,1,l,1);
        trA_L(l/10-1,t) = tr_L;

        %no Nystrom, 2 Hutch l/2 Lanczos
        [~,tr_HL] = EST(A,mi,2,l/2,1);
        trA_HL(l/10-1,t) = tr_HL;

        %no Nystrom, 5 Hutch l/5 Lanczos
        [~,tr_5HL] = EST(A,mi,5,l/5,1);
        trA_5HL(l/10-1,t) = tr_5HL;
        
    end

end

%means
mtr1 = 1/T * sum(tr1');
mtr2 = 1/T * sum(tr2');
mtr5 = 1/T * sum(tr5');
mitr1 = 1/T * sum(itr1');
mitr2 = 1/T * sum(itr2');
mitr5 = 1/T * sum(itr5');
mtrA_L = 1/T * sum(trA_L');
mtrA_HL = 1/T * sum(trA_HL');
mtrA_5HL = 1/T * sum(trA_5HL');

figure(2)
plot(mvecs,abs(mitr5-trA)./abs(trA),'k')
hold on
plot(mvecs,abs(mitr2-trA)./abs(trA),'m')
hold on
plot(mvecs,abs(mtr5-trA)./abs(trA),'r')
hold on
plot(mvecs,abs(mtr1-trA)./abs(trA),'g')
hold on
plot(mvecs,abs(mtr2-trA)./abs(trA),'y')
hold on
plot(mvecs,abs(mtrA_L-trA)./abs(trA),'c')
hold on
plot(mvecs,abs(mtrA_HL-trA)./abs(trA),'b')

legend('inv Nystrom 5Hutch 2Lanczos', 'inv Nystrom 2Hutch 5Lanczos','Nystrom 5Hutch 2Lanczos', 'Nystrom 1Hutch 5Lanczos', 'Nystrom 2Hutch 5Lanczos', '1Hutch l-Lanczos', '2Hutch 0.5l-Lanczos')
xlabel('matvecs')
ylabel('error')
title('relative error wrt matvecs')

% % intervals
% Max_tr5 = max(tr5');
% Min_tr5 = min(tr5');
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
% plot(mvecs,Max_tr1,'g')
% hold on
% GG = plot(mvecs,mtr1,'o-g')
% hold on
% plot(mvecs,Min_tr1,'g')
% hold on
% plot(mvecs,Max_tr2,'y')
% hold on
% YY = plot(mvecs,mtr2,'o-y')
% hold on
% plot(mvecs,Min_tr2,'y')
% hold on
% plot(mvecs,Max_trA_L,'c')
% hold on
% CC = plot(mvecs,mtrA_L,'o-c')
% hold on
% plot(mvecs,Min_trA_L,'c')
% hold on
% plot(mvecs,Max_trA_HL,'b')
% hold on
% BB = plot(mvecs,mtrA_HL,'o-b')
% hold on
% plot(mvecs,Min_trA_HL,'b')
% 
% legend([RR, GG, YY, CC, BB],'Nystrom 5Hutch 2Lanczos', 'Nystrom 1Hutch 5Lanczos', 'Nystrom 2Hutch 5Lanczos', '1Hutch l-Lanczos', '2Hutch 0.5l-Lanczos')
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