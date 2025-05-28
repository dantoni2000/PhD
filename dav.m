% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off

T = 50;

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

    case 1
        n = 100;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        % mi = 1;
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
        n = 100;
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
        n = 100;
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
        n = 5000; mi = 1;
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
        A = zeros(n,n);
        l=1; h=50;
        for j=1:40
            x = sprand(n,1,0.01);
            A = A + (h/j^2*x) * x';
        end
        for j=41:300
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
end

trA = sum(log(diag(G+mi*eye(n,n))),"all");

for t=1:T

    l = 100;
    [UBig,LhatBig] = Nystrom(A,l);
    PBig = (mi)^0.5*UBig*(LhatBig+mi*eye(l))^-0.5*UBig' + (eye(n) - UBig*UBig');
    trPBig = sum(log(diag(LhatBig+mi*eye(l,l))),"all");
    
    l1 = 50;
    [U1,Lhat1] = Nystrom(A,l1);  ll1 = 0*Lhat1(l1,l1);
    P1 = (mi)^0.5*U1*(Lhat1+mi*eye(l1))^-0.5*U1' + (eye(n) - U1*U1');
    trP1= sum(log(1/(ll1+mi)*diag(Lhat1+mi*eye(l1,l1))),"all");

    l2 = 10;
    [U2,Lhat2] = Nystrom(A,l2);  ll2 = 0*Lhat2(l2,l2);
    P2 = (mi)^0.5*U2*(Lhat2+mi*eye(l2))^-0.5*U2' + (eye(n) - U2*U2');
    trP2 = sum(log(1/(ll2+mi)*diag(Lhat2+mi*eye(l2,l2))),"all");

    for s=20:20:200

        mvecs(s/20) = 10*s + 100;

        % 2s samples 5 its Lanczos
        [~,PtrBig] = Preconditioned_HUTCH(A,mi,PBig,2*s,5,1,1);
        trBig(s/20,t) = trPBig+ PtrBig;

        % [~,PtrBig2] = Preconditioned_HUTCH(A,mi,PBig,s,10,1,1);
        % trBig2(s/20,t) = trPBig+ PtrBig2;       
    
        % 2s + 4 samples 5 its Lanczos
        [~,Ptr1] = Preconditioned_HUTCH(A,mi,P1,2*s + 10,5,1,1);
        tr1(s/20,t) = trP1 + Ptr1;

        % [~,Ptr12] = Preconditioned_HUTCH(A,mi,P1,s + 5,10,1,1);
        % tr12(s/20,t) = trP1 + Ptr12;     

        % 2s + 8 samples 5 its Lanczos
        [~,Ptr2] = Preconditioned_HUTCH(A,mi,P2,2*s + 18,5,1,1);
        tr2(s/20,t) = trP2 + Ptr2;

        % [~,Ptr22] = Preconditioned_HUTCH(A,mi,P2,s + 9,10,1,1);
        % tr22(s/20,t) = trP2 + Ptr22;
        
    end

end

%means
mtrBig = 1/T * sum(trBig');
% mtrBig2 = 1/T * sum(trBig2');
mtr1 = 1/T * sum(tr1');
% mtr12 = 1/T * sum(tr12');
mtr2 = 1/T * sum(tr2');
% mtr22 = 1/T * sum(tr22');

figure(2)
plot(mvecs,abs(mtrBig-trA)./abs(trA),'m')
hold on
% plot(mvecs,abs(mtrBig2-trA)./abs(trA),'r')
% hold on
plot(mvecs,abs(mtr1-trA)./abs(trA),'g')
hold on
% plot(mvecs,abs(mtr12-trA)./abs(trA),'y')
% hold on
plot(mvecs,abs(mtr2-trA)./abs(trA),'c')
hold on
% plot(mvecs,abs(mtr22-trA)./abs(trA),'b')

legend('Big Nystrom', 'Mid Nystrom', 'Small Nystrom')
xlabel('matvecs')
ylabel('error')
title('relative error wrt matvecs')