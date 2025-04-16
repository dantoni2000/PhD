% 09 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off

T = 5; dl = 10; S = 100;

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
        n = 2000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-3;
        G = diag(1./(g).^2);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 2
        n = 2000;
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
        n = 2000;
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
        n = 2000;
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
        n = 2000;
        mi = 1e-4;
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

for t=1:T
    
    flopsHutch = zeros(dl,1); flopsPreHutch = zeros(dl,1); 

    for s=10:10:S
        
        % Nystrom su A
        [U,Lhat] = Nystrom(A,s);  ll = Lhat(s,s);
        % P05 = sqrtm((ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U') );
        % Pnorm = U*(Lhat+mi*eye(l))^-0.5*U' + (ll+mi)^-0.5*(eye(n) - U*U');
        P = (ll+mi)^0.5*U*(Lhat+mi*eye(s))^-0.5*U' + (eye(n) - U*U');        
            
        flopsPreHutch(s/10,1) = flopsPreHutch(s/10,1) + s^2*(s+n);
        % trPREC(l/10,1) = trace(logm((P*P)^-1));
        trPREC(s/10,1) = log((ll+mi)^(-s)*prod(diag(Lhat+mi*eye(s))));
    
        flopsPreHutch(s/10,1) = flopsPreHutch(s/10,1) + 3*s;
        % true_value(l/10,1) = trace(logm(P05*(A+mi*eye(n))*P05));

        %Hutchinson sulla matrice precondizionata
        [Pits(s/10,t),Ptr] = Preconditioned_HUTCH(A,mi,P,s,500);
        flopsPreHutch(s/10,1) = flopsPreHutch(s/10,1) + s*(Pits(s/10)*n^2);
        PreHutch(s/10,1,t) = trPREC(s/10,1) + Ptr;
        %stimatore direttamente su A
                
        [its(s/10,t),tr] = EST(A,mi,s,500);
        HutchA(s/10,t) = tr;
        flopsHutch(s/10,1) = flopsHutch(s/10,1) + s*(its(s/10)*n^2);
        %stimatore deterministico
        brutal_estimate(s/10,t) = n*log(ll); %+G(n,n)

        %stime finali
        BrutEst(s/10) = trPREC(s/10) + brutal_estimate(s/10);

    end

        figure(2)
        plot(10:10:S,PreHutch(:,t),'r')
        hold on
        plot(10:10:S,HutchA(:,t),'o-g')
        hold on
        plot(10:10:S,trA*ones(S/10,1),'c')
        legend('PreHutch++','est(log(A+$\mu$ I))','trlog(A+$\mu$ I)','interpreter','latex')
        xlabel('samples for the Nystrom preconditioner')
        ylabel('trace computation')

        figure(3)
        plot(10:10:S,flopsHutch(:,1),'o-g')
        hold on
        plot(10:10:S,flopsPreHutch(:,1),'r')
        legend('FLOPS est(log(A+$\mu$ I))','FLOPS PreHutch++','interpreter','latex')
        xlabel('samples for the Nystrom preconditioner')
        ylabel('flops')
end