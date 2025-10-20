clear all
close all
warning off

T = 10;
S = [];

decadimento=3;

switch decadimento

    case 0
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1;
        G = 10^2*diag((-g+g(n,1)+1)./g(n,1));
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
        G = 10^2*diag(1./(g).^.85);
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
        G = 10^7*diag(exp(-.125*g));
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
        G=10^8*diag(diag(G));
        % G = G + .1*eye(n);
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
        G = 10^3*sparse(diag(exp(-(1:n)/100)));
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
        % A = A + .5*eye(n);
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
        A = 10^3*A;
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
nA = sum(log(1+diag(G)).^2).^0.5;
nTr = zeros(21,1); stDevnTr  = zeros(21,1);
nFro = zeros(21,1); stDevnFro  = zeros(21,1);
nNoPrec = zeros(21,1); stDevnNoPrec  = zeros(21,1);

dl = 20; s = 10;

    for l=10:10:10*dl
        
        % Nystrom su A
        [U,Lhat] = Nystrom(A,l);  ll = Lhat(l,l);
        % P05 = sqrtm((ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U') );
        % Pnorm = U*(Lhat+mi*eye(l))^-0.5*U' + (ll+mi)^-0.5*(eye(n) - U*U');
        P = (ll+mi)^0.5*U*(Lhat+mi*eye(l))^-0.5*U' + (eye(n) - U*U');
        % plot(eigs(P*(A+mi*eye(n))*P));
        trPREC(l/10,1) = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l))));
        % trPREC(l/10,1) = trace(logm((P*P)^-1));

        %Hutchinson sulla matrice precondizionata
        figure(1+l/10)
        [Pits(l/10),Ptr,h] = Preconditioned_HUTCH(A,mi,P,s,500,0,0);
        PreHutch(l/10,s/10) = trPREC(l/10,1) + Ptr;
        cPAP = 2*sqrt(cond(P*A*P)+1)*log(2*cond(P*A*P));
        rhoPAP = ((sqrt(cond(P*A*P)+1)-1)/(sqrt(cond(P*A*P)+1)+1))^2;
        % bound teorico
        % b = semilogy(cPAP*rhoPAP.^[1:fix(Pits(l/10))+1],'b');
        
        hold on
        [its(l/10),tr,k] = EST(A,mi,s,500,0);
        HutchA(l/10,s/10) = tr;
        % f = semilogy(cA*rhoA.^[1:fix(its(l/10))+1],'m');

        legend([h k],'$\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v - v^T \log(T^{(P)}_m) v \right\|$','$\left\| v^T \log(A + \mu I) v - v^T \log(T_m) v \right\|$','interpreter','latex','FontSize',12)
        %legend([h k],'$\frac{\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v - v^T \log(T^{(P)}_m) v \right\|}{\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v \right\|}$','$\frac{\left\| v^T \log(A + \mu I) v - v^T \log(T_m) v \right\|}{\left\| v^T \log(A + \mu I) v \right\|}$','interpreter','latex','FontSize',12)
        xlabel('Lanczos Iterations')
        ylabel('Error')
    end