% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off

T = 1; dl = 20; S = 10;
MatVecsEST = zeros(dl,1); MatVecsHutch = zeros(dl,1); 

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
[its0,tr0] = EST(A,mi,S,500,1);
MatVecsEST0 = its0*S;

for t=1:T
    for l=10:10:10*dl
    
        % Nystrom su A
        [U,Lhat] = Nystrom(A,l);  ll = Lhat(l,l);
        % P05 = sqrtm((ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U') );
        P = (ll+mi)^0.5*U*(Lhat+mi*eye(l))^-0.5*U' + (eye(n) - U*U');
        
        MatVecsHutch(l/10) = MatVecsHutch(l/10) + l;
        % trPREC(l/10,1) = trace(logm((P*P)^-1));
        trPREC(l/10,1) = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l))));
        % true_value(l/10,1) = trace(logm(P05*(A+mi*eye(n))*P05));
        
        %Hutchinson sulla matrice precondizionata
        [Pits(l/10),Ptr] = Preconditioned_HUTCH(A,mi,P,S,500,1,1);
        MatVecsHutch(l/10) = MatVecsHutch(l/10) + Pits(l/10)*S;
        PreHutch(l/10,1) = trPREC(l/10,1) + Ptr;
        %stimatore direttamente su A
        Z = fix(MatVecsHutch(l/10)/fix(its0))+1;
        [its(l/10),tr] = EST(A,mi,Z,500,1);
        MatVecsEST(l/10) = its(l/10)*Z;
        HutchA(l/10,1) = tr;


        %stimatore deterministico
        brutal_estimate(l/10) = n*log(ll); %+G(n,n)

        %stime finali
        BrutEst(l/10) = trPREC(l/10) + brutal_estimate(l/10);        
    end
    [sMatVecsEST,ord] = sort(MatVecsEST);


    figure(2)
    plot(sMatVecsEST,PreHutch(ord)-trA,'r')
    hold on
    plot(sMatVecsEST,HutchA(ord)-trA,'o-g')
    legend('PreHutch++','est(log(A+$\mu$ I))','trlog(A+$\mu$ I)','interpreter','latex')

end