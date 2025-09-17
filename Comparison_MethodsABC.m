clear all
close all
warning off

T = 10;

decadimento=5;

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
        G = 10^3*diag(1./sqrt(g));
        A = Q*G*Q';
        % A = G;

    case 3
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10^3*diag(exp(-.25*g));
        % G = 1000*diag(exp(-.1*g));
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
        G=10^6*diag(diag(G));
        A = Q*G*Q'; 
        % A = G;
        
    case 6
        n = 1000; mi = 1;
        G = 10^5*(diag((1:n).^(-2)));
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
        g = linspace(1,n,n)';
        mi = 1;
        G = eye(n);
        A = Q*G*Q';
        
end

mi = 1;
trA = sum(log(diag(G+mi*eye(n,n))),"all");
nA = norm(A,'fro');
nTr = zeros(19,1); stDevnTr  = zeros(19,1);
nFro = zeros(19,1); stDevnFro  = zeros(19,1);
nNoPrec = zeros(19,1); stDevnNoPrec  = zeros(19,1);
nxNys = zeros(19,1);

for l=15:10:195
        
    mvecs((l-5)/10 ) = l;
        
    for t=1:T
        % Nystrom grande su A
        [UBig,LhatBig] = PinvNystrom(A,l);
        nTr((l-5)/10) = nTr((l-5)/10) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA);
        stDevnTr((l-5)/10) = stDevnTr((l-5)/10) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA).^2;

        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-5;
        [~,trr] = Nystrom_HUTCH(A,mi,l1,1,5,1,1);
        nFro((l-5)/10) = nFro((l-5)/10) + abs(trr - trA);
        
        %Nystrom con n Hutch
        [~,trrNoPrec] = EST(A,mi,5,l/5,1);
        nNoPrec((l-5)/10) = nNoPrec((l-5)/10) + abs(trrNoPrec - trA);

        % % Per ora xnystrace con f esatta. Come si mette f in xnystrace?
        % [xnys,~] = xnystrace(logm(A+eye(n)), l1);
        % nxNys((l-5)/10) = nxNys((l-5)/10) + abs(xnys - trA);

    end

end

nFro = 1/T * nFro;
nTr = 1/T * nTr;
nNoPrec = 1/T * nNoPrec;
nxNys = 1/T * nxNys;

for tMV = 10:10:190

    for s = 4:2:tMV-4
        k = s;
        p = tMV - s;
        
        % boundTr(:,s/2-1) = (1 + k/(p-1)).* norm(G(k+1:n,k+1:n).^0.5 ,'fro')^2;
        % boundFro(:,s/2-1) = (1 + sqrt(k/(p-1))).^2.* norm(G(k+1:n,k+1:n),'fro') ;
        % boundFro(:,s/2-1) = norm(G(k+1:n,k+1:n),'fro') + k/(p-1) .* norm(G(k+1:n,k+1:n).^0.5 ,'fro')^2 ;
        % boundSpec(:,s/2-1) = (1 + 2*k/(p-1)).* norm(G(k+1:n,k+1:n)) + (2*exp(1)^2*(k+p)/(p^2 - 1)) .* norm(G(k+1:n,k+1:n).^0.5 ,'fro')^2;
        % conjFro(:,s/2-1) = (1 + k/(p-1)).* norm(G(k+1:n,k+1:n),'fro');

        squareboundTr(:,s/2-1) = sqrt( 2* (1+sqrt(exp(1)^5/2)*((k+5)/(p+1))^2 )* norm(G(k+5+1:n,k+5+1:n).^0.5,'fro')^4 + 2 * (k+5)*(k+5+p-1)/(p*(p-1)*(p-3)) * norm(G(k+5+1:n,k+5+1:n),'fro')^2);
        squareboundFro(:,s/2-1) = sqrt( 4* (1+sqrt(exp(1)^5/2)*(k/(p+1))^2 + k*(k+p-1)/(p*(p-1)*(p-3)))* norm(G(k+1:n,k+1:n),'fro')^2 + k*(k+p-1)/(p*(p-1)*(p-3)) * norm(G(k+1:n,k+1:n).^0.5,'fro')^4);
        % squareboundSpec(:,s/2-1) = sqrt( 2* (1+12*sqrt(exp(1)^5/2)*(k/(p+1))^2)* norm(G(k+1:n,k+1:n))^2 + (12*exp(1)^4)*(k+p).^2/((p+1)^3*(p-3)) * norm(G(k+1:n,k+1:n).^0.5,'fro')^4);
    end

    BestSqrtTr(tMV/10) = 1/mi * min(squareboundTr');
    BestSqrtFro(tMV/10) = 1/mi * min(squareboundFro');
    BoundNoPrec(tMV/10) = 1/mi * sqrt(10/tMV)*nA;
    % BestSqrtSpec(tMV/10) = min(squareboundSpec');

    BestRkTr(tMV/10) = 1/mi * norm(diag(G(tMV+1:n,tMV+1:n)).^0.5,'fro').^2;
    BestRkFro(tMV/10) = 1/mi * norm(diag(G(tMV+1:n,tMV+1:n)),'fro');
    BestRkSpec(tMV/10) = 1/mi * norm((G(tMV+1:n,tMV+1:n)));

    check(tMV/10) = sum(diag(G(k+1:k+5,k+1:k+5))) / sum(diag(G(k+5+1:n,k+5+1:n)));

end
% figure(100)
% plot(diag(A))

% tipo qualcosa che dipende dal gap ex 1/(1+gap*l) (1+l) ||A-A_k+p+l||
% for jj = 1:20
%     ratio(jj) = nFro(end)./norm(diag(G(tMV+1+jj:n,tMV+1+jj:n)),'fro');
% end
% figure(4)
% plot(ratio)

figure(4)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(diag(G))
xlabel('$n$','interpreter','Latex')
ylabel('eigenvalues')
title('Eigenvalues of the matrix')
legend('$\lambda(A)$','interpreter','Latex')

figure(4)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(mvecs,nTr,'-b')
hold on
semilogy(mvecs,BestSqrtTr','-*c')
hold on
semilogy(mvecs,nFro, '-r');
hold on
semilogy(mvecs,BestSqrtFro','-*m')
hold on
semilogy(mvecs,nNoPrec, 'Color', [0.6350 0.0780 0.1840], 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 8)
hold on
semilogy(mvecs,BoundNoPrec','-*', 'Color', [0.9290 0.6940 0.1250], 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 8)
hold on
% semilogy(mvecs,nxNys, '-k')
xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for the two strategies')
legend('error (A)', 'bound (A)', 'error (B)','bound (B)', 'error (C)', 'bound (C)','fontsize',18)