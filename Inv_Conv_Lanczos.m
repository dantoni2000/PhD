% 08 04 2025 convergence of  Lanczos for the original and for the
% preconditioned matrix

clear all
close all
warning off

T = 1; dl = 20; s = 10;

decadimento=3;

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
        mi = 5e-6;
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
        mi = 1;
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
end

trA = sum(log(diag(G+mi*eye(n,n))),"all");
kA = cond(A);
cA = 2*sqrt(kA+1)*log(2*kA);
rhoA = ((sqrt(kA+1)-1)/(sqrt(kA+1)+1))^2;

for t=1:T 
        
    for l=10:10:10*dl
        
        % Nystrom su inv(A)
        [iU,iLhat] = Nystrom(inv(A+1e-14*eye(n)),l);  iLhat=inv(iLhat); ill = iLhat(l,l);
        iP = (ill+mi)^0.5*iU*(iLhat+mi*eye(l))^-0.5*iU' + (eye(n) - iU*iU');
        % trPRECi(l/10,1) = log((ill+mi)^(-l)*prod(diag(iLhat+mi*eye(l,l))));

        % Nystrom su A
        [U,Lhat] = Nystrom(A,l);  ll = 0 * Lhat(l,l);
        P = (ll+mi)^0.5*U*(Lhat+mi*eye(l))^-0.5*U' + (eye(n) - U*U');
        % plot(eigs(P*(A+mi*eye(n))*P));
        % trPREC(l/10,1) = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l))));

        % combinazione dei due Nystrom
        V = [U(:,1:l/2) iU(:,1:l/2)];
        Lambdahat = diag([diag(Lhat(1:l/2,1:l/2)); diag(iLhat(1:l/2,1:l/2))]);
        Q = mi^0.5*V*(Lambdahat+mi*eye(l))^-0.5*V' + (eye(n) - V*V');

        figure(1+l/10)

        %Hutchinson sulla matrice precondizionata con Nystrom sull'inversa
        [iPits(l/10),iPtr,ih] = Preconditioned_HUTCH(A,mi,iP,s,500,1,0);
        % PreHutch(l/10,s/10) = trPREC(l/10,1) + iPtr;
        ciPAP = 2*sqrt(cond(iP*(A+mi*eye)*iP)+1)*log(2*cond(iP*(A+mi*eye)*iP));
        rhoiPAP = ((sqrt(cond(iP*(A+mi*eye)*iP)+1)-1)/(sqrt(cond(iP*(A+mi*eye)*iP)+1)+1))^2;
        % bound teorico
        % b = semilogy(cPAP*rhoPAP.^[1:fix(Pits(l/10))+1],'b');        
        hold on

        %Hutchinson sulla matrice precondizionata con Nystrom
        [Pits(l/10),Ptr,h] = Preconditioned_HUTCH(A,mi,P,s,500,0,0);
        cPAP = 2*sqrt(cond(P*(A+mi*eye)*P)+1)*log(2*cond(P*(A+mi*eye)*P));
        rhoPAP = ((sqrt(cond(P*(A+mi*eye)*P)+1)-1)/(sqrt(cond(P*(A+mi*eye)*P)+1)+1))^2;
        hold on

        %Hutchinson sulla matrice precondizionata con media dei due Nystrom
        [Qits(l/10),Qtr,md] = Preconditioned_HUTCH(A,mi,Q,s,500,2,0);
        cQAQ = 2*sqrt(cond(Q*(A+mi*eye)*Q)+1)*log(2*cond(Q*(A+mi*eye)*Q));
        rhoQAQ = ((sqrt(cond(Q*(A+mi*eye)*Q)+1)-1)/(sqrt(cond(Q*(A+mi*eye)*Q)+1)+1))^2;
        hold on

        [its(l/10),tr,k] = EST(A,mi,s,500,0);
        % HutchA(l/10,s/10) = tr;
        % f = semilogy(cA*rhoA.^[1:fix(its(l/10))+1],'m');

        legend([ih h md k],'$\left\| v^T \log(P_i^{\frac{1}{2}} (A + \mu I) P_i^{\frac{1}{2}}) v - v^T \log(T^{(P_i)}_m) v \right\|$','$\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v - v^T \log(T^{(P)}_m) v \right\|$','$\left\| v^T \log(Q^{\frac{1}{2}} (A + \mu I) Q^{\frac{1}{2}}) v - v^T \log(T^{(Q)}_m) v \right\|$','$\left\| v^T \log(A + \mu I) v - v^T \log(T_m) v \right\|$','interpreter','latex','FontSize',12)
        %legend([h k],'$\frac{\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v - v^T \log(T^{(P)}_m) v \right\|}{\left\| v^T \log(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}}) v \right\|}$','$\frac{\left\| v^T \log(A + \mu I) v - v^T \log(T_m) v \right\|}{\left\| v^T \log(A + \mu I) v \right\|}$','interpreter','latex','FontSize',12)
        xlabel('Lanczos Iterations')
        ylabel('Error')
        title('convergence of Lanczos method for different preconditioners')
    end

    figure
    plot(iPits,'c')
    hold on
    plot(Pits,'r')
    hold on
    plot(Qits,'b')
    hold on
    plot(its,'g')
    legend('$(P_i^{\frac{1}{2}} (A + \mu I) P_i^{\frac{1}{2}})$','$(P^{\frac{1}{2}} (A + \mu I) P^{\frac{1}{2}})$','$(Q^{\frac{1}{2}} (A + \mu I) Q^{\frac{1}{2}})$','$(A + \mu I)$','interpreter','latex','FontSize',12)
    xlabel('rank of preconditioner')
    ylabel('iterations')
    title('number of iteration of the Lanczos method for different preconditioners wrt their rank')
end