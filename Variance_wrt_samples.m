% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off

T = 1; dl = 20; S = 100;
flopsHutch = zeros(dl,S/10); flopsPreHutch = zeros(dl,S/10); 

decadimento=0;

switch decadimento
    case 0
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 0;
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
        mi = 1e-3;
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
MeanTotVarHutch = zeros(dl,1); MeanTotVarPreHutch = zeros(dl,1);

for t=1:T
    for l=10:10:10*dl
    
        % Nystrom su A
        [U,Lhat] = Nystrom(A,l);  ll = Lhat(l,l);
        % P05 = sqrtm((ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U') );
        P = (ll+mi)^0.5*U*(Lhat+mi*eye(l))^-0.5*U' + (eye(n) - U*U');
        
        flopsPreHutch(l/10,1:S/10) = flopsPreHutch(l/10,1) + l^2*(l+n);
        % trPREC(l/10,1) = trace(logm((P*P)^-1));
        trPREC(l/10,1) = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l))));

        flopsPreHutch(l/10,1:S/10) = flopsPreHutch(l/10,1) + 3*l;
        % true_value(l/10,1) = trace(logm(P05*(A+mi*eye(n))*P05));

        for s=10:10:S
            %Hutchinson sulla matrice precondizionata
            [Pits(s/10),Ptr] = Preconditioned_HUTCH(A,mi,P,s,500);
            flopsPreHutch(l/10,s/10) = flopsPreHutch(l/10,s/10) + s*(Pits(s/10)*n^2);
            PreHutch(l/10,s/10) = trPREC(l/10,1) + Ptr;
            %stimatore direttamente su A
            [its(s/10),tr] = EST(A,mi,s,500);
            HutchA(l/10,s/10) = tr;
            flopsHutch(l/10,s/10) = flopsHutch(l/10,s/10) + s*(its(s/10)*n^2);
        end
        
        %total variation of estimator and PreHutch
        tvh(l/10,:) = HutchA(l/10,2:end)-HutchA(l/10,1:end-1);
        TotVarHutch(l/10,t) = sum(abs(tvh(l/10,:)));
        MeanTotVarHutch(l/10,1) = MeanTotVarHutch(l/10) + TotVarHutch(l/10,t);
        tvph(l/10,:) = PreHutch(l/10,2:end)-PreHutch(l/10,1:end-1);
        TotVarPreHutch(l/10,t) = sum(abs(tvph(l/10,:)));
        MeanTotVarPreHutch(l/10,1) = MeanTotVarPreHutch(l/10) + TotVarPreHutch(l/10,t);

        %stimatore deterministico
        brutal_estimate(l/10) = n*log(ll); %+G(n,n)

        %stime finali
        BrutEst(l/10) = trPREC(l/10) + brutal_estimate(l/10);

        figure(1+l/10)
        plot(10:10:S,PreHutch(l/10,:),'r')
        hold on
        plot(10:10:S,HutchA(l/10,:),'o-g')
        hold on
        plot(10:10:S,trA*ones(S/10,1),'c')
        legend('PreHutch++','est(log(A+$\mu$ I))','trlog(A+$\mu$ I)','interpreter','latex')
        
        figure(100+l/10)
        plot(10:10:S,flopsHutch(l/10,:),'o-g')
        hold on
        plot(10:10:S,flopsPreHutch(l/10,:),'r')
        legend('FLOPS est(log(A+$\mu$ I))','FLOPS PreHutch++','interpreter','latex')
    end
    
    % figure(3)
    % plot([10:10:10*dl],PreHutchEst/trA,'r')
    % hold on
    % plot([10:10:10*dl],PreHutchEstInv/trA,'-p')
    % hold on
    % plot([10:10:10*dl],HutchEst/trA,'b')
    % hold on
    % plot([10:10:10*dl],HutchEstInv/trA,'-v')
    % hold on
    % plot([10:10:10*dl],BrutEst/trA,'g')
    % hold on
    % plot([10:10:10*dl],BrutEstInv/trA,'y')
    % hold on
    % plot([10:10:10*dl],HutchA/trA,'o-')
    % hold on
    % ylim([0.8,1.2])
    % legend('ratio PreHutch++','ratio PreHutchInv++','ratio Hutch++','ratio HutchInv++','ratio BrutalEst','ratio BrutalEstInv','ratio Hutch(A)')
end

figure(400*T)
plot(10:10:10*dl,1/T*MeanTotVarHutch,'g')
hold on
plot(10:10:10*dl,1/T*MeanTotVarPreHutch,'r')
xlabel('rank of Nystrom preconditioner')
ylabel('mean of total variance')
legend('Total Variance Estimate','Total Variance PreHutch++','interpreter','latex')