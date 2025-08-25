clear all
close all
warning off

T = 30;

decadimento=8;

switch decadimento

    case 0
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10^3*diag((-g+g(n,1)+1)./g(n,1));
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
        G = 10^3*diag(1./(g).^.85);
        A = Q*G*Q';
        % A = G;

    case 2
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10^3 * diag(1./sqrt(g));
        A = Q*G*Q';
        % A = G;

    case 3
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10^4*diag(exp(-.25*g));
        A = Q*G*Q';
        % A = G;

    case 4
        n = 1000;
        mi = 1;
        kernel = @(x,y) ((norm(x)/n)^10 + (norm(y)/n)^10)^(0.1);
        data_matrix = orth(randn(n,n))*diag(1:n);
        A = build_kernel_matrix(data_matrix,kernel);
        [Q,G] = svd(A);
        G=10^6*diag(diag(G));
        A = Q*G*Q'; 
        % A = G;
        
    case 5
        n = 1000;
        mi = 1;
        alpha = 1; nu = 3/2; % nu = 13/2;
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
        G=10^6*diag(diag(G));
        A = Q*G*Q'; 
        % A = G;
        
    case 6
        n = 1000; mi = 1;
        G = 10^5*sparse(diag((1:n).^(-2)));
        % Q = gallery('orthog',n,1);
        [Q,~] = qr(randn(n,n));
        A = Q*G*Q';  
        % A = G;
        
    case 7
        n = 1000; mi = 1;
        G = 10^4*sparse(diag(exp(-(1:n)/100)));
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
        A = 10^6 * A;
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
        A = 10^3 * A;
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
nFro2 = zeros(9,1);
nFro3 = zeros(9,1);
nFro4 = zeros(9,1);
nNoPrec = zeros(9,1);

for l=25:10:105
        
    mvecs((l-15)/10 ) = l;
        
    for t=1:T
               
        %Nystrom con 1 Hutch 5 Lanczos
        [UBig,LhatBig] = PinvNystrom(A,l);
        nTr((l-15)/10) = nTr((l-15)/10) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA);


        l1 = l-5;
        [~,trr] = Nystrom_HUTCH(A,mi,l1,1,5,1,1);
        nFro((l-15)/10) = nFro((l-15)/10) + abs(trr - trA);

        l2 = l-10;
        [~,trr2] = Nystrom_HUTCH(A,mi,l2,2,5,1,1);
        nFro2((l-15)/10) = nFro2((l-15)/10) + abs(trr2 - trA);
        
        l3 = l-15;
        [~,trr3] = Nystrom_HUTCH(A,mi,l3,3,5,1,1);
        nFro3((l-15)/10) = nFro3((l-15)/10) + abs(trr3 - trA);

        l4 = l-20;
        [~,trr4] = Nystrom_HUTCH(A,mi,l4,4,5,1,1);
        nFro4((l-15)/10) = nFro4((l-15)/10) + abs(trr4 - trA);

        %Nystrom con n Hutch
        [~,trrNoPrec] = EST(A,mi,l/5,5,1);
        nNoPrec((l-15)/10) = nNoPrec((l-15)/10) + abs(trrNoPrec - trA);

    end

end

nTr = 1/T * nTr;
nFro = 1/T * nFro;
nFro2 = 1/T * nFro2;
nFro3 = 1/T * nFro3;
nFro4 = 1/T * nFro4;
nNoPrec = 1/T * nNoPrec;

figure(3)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(3)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nTr,'m')
hold on
semilogy(mvecs,nFro,'r')
hold on
semilogy(mvecs,nFro2,'b')
hold on
semilogy(mvecs,nFro3,'g')
hold on
semilogy(mvecs,nFro4,'c')
hold on
semilogy(mvecs,nNoPrec,'k')
% hold on
% semilogy(mvecs,BestRkTr','-*b')
% hold on
% semilogy(mvecs,BestLowerFro','-dr')
% hold on
% semilogy(mvecs,BestLowerTr','-db')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for logtrace')
legend('0 sample', '1 sample', '2 sample', '3 sample','4 sample', 'tutti sample')