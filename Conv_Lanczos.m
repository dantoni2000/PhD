% 08 04 2025 convergence of  Lanczos for the original and for the
% preconditioned matrix

clear all
close all
warning off

T = 1; dl = 20; s = 10;

decadimento=5;

switch decadimento
    case 0
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-2;
        G = diag((-g+g(n,1)+1)./g(n,1));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 1
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-4;
        G = diag(1./(g).^2);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on
        semilogy(diag(G)+mi*ones(n))

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

    case 3
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-3;
        G = diag(exp(-g));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

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

    case 5
        n = 1000;
        mi = 5e-6;
        alpha = 1; nu = 5/2;
        kernel = @(x,y) sqrt(pi)*((alpha*norm(x-y))^(nu)*besselk(nu,alpha*norm(x-y)))/(2^(nu-1)*alpha^(2*nu)*gamma(nu+0.5));
        data_matrix = 1/n*randn(n,n);

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
end

trA = sum(log(diag(G+mi*eye(n,n))),"all");
kA = cond(A);
cA = 2*sqrt(kA+1)*log(2*kA);
rhoA = ((sqrt(kA+1)-1)/(sqrt(kA+1)+1))^2;

for t=1:T 
        
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

end