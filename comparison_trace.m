% 30 03 2025 Test Hutchinson: in quali ipotesi lo stimatore gratuito funziona meglio
% DATI SINTETICI

clear all
close all
warning off

T = 1; dl = 10;

decadimento=3;

switch decadimento
    case 0
        n = 500;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 0;
        G = diag((-g+g(n,1)+1)./g(n,1));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        
end

switch decadimento
    case 1
        n = 500;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 0;
        G = diag(1./(g).^2);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        
end

switch decadimento
    case 2
        n = 500;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 0;
        G = diag(1./sqrt(g));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
end

switch decadimento
    case 3
        n = 500;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 0;
        eta = 1e-6;
        G = diag(exp(-g)+eta*ones(n,1));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
end

switch decadimento
    case 4
        n = 500;
        mi = 0;
        eta = 1e-6;
        kernel = @(x,y) ((norm(x)/n)^10 + (norm(y)/n)^10)^(0.1);
        data_matrix = orth(randn(n,n))*diag(1:n);
        A = build_kernel_matrix(data_matrix,kernel);
        [Q,G] = svd(A);
        G=diag(diag(G)+eta*ones(n,1));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
end

trA = sum(log(diag(G+mi*eye(n,n))),"all");

for t=1:T
    for l=10:10:10*dl
    
        % Nystrom su A
        [U,Lhat] = Nystrom(A,l);  ll = Lhat(l,l); %ll = Lhat(fix(l/2)+1,fix(l/2)+1);
        P05 = sqrtm((ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U') );
        P = (ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U');
        
        trPREC(l/10) = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l))));

        true_value(l/10) = trace(logm(P05*(A+mi*eye(n))*P05));

        %Nystrom su inv(A)
        [Uinv,Lhatinv] = Nystrom(A^-1,l); Lhatinv = Lhatinv^-1; llinv = Lhatinv(l,l); %llinv = Lhatinv(fix(l/2)+1,fix(l/2)+1); 
        B05 = sqrtm((llinv+mi)*Uinv*(Lhatinv+mi*eye(l))^-1*Uinv' + (eye(n) - Uinv*Uinv') );
        B = (llinv+mi)*Uinv*(Lhatinv+mi*eye(l))^-1*Uinv' + (eye(n) - Uinv*Uinv');
        
        trPRECinv(l/10) = log((llinv+mi)^(-l)*prod(diag(Lhatinv+mi*eye(l,l))));
    
        true_valueinv(l/10) = trace(logm(B05*(A+mi*eye(n))*B05));
    
        %Hutchinson ''base''
        Hutch(l/10) = 0;
        Hutchinv(l/10) = 0;
        
        for j=1:l
            
            x=randn(n,1);
            Hutch(l/10) = Hutch(l/10) + Lanczos_log((P05*(A+mi*eye(n))*P05),x,5);
            Hutchinv(l/10) = Hutchinv(l/10) + Lanczos_log((B05*(A+mi*eye(n))*B05),x,5);
            
        end

        Hutch(l/10)= 1/l*Hutch(l/10);
        Hutchinv(l/10)= 1/l*Hutchinv(l/10);

        %Hutchinson sulla matrice precondizionata
        PreHutch(l/10) = Preconditioned_HUTCH(A,mi,P,l,5);
        PreHutchinv(l/10) = Preconditioned_HUTCH(A,mi,B,l,5);

        %stimatore deterministico
        brutal_estimate(l/10) = n*log(ll); %+G(n,n)
        brutal_estimate_inv(l/10) = n*log(llinv); %+G(n,n)
        
        %Hutchinson direttamente su A
        HutchA(l/10) = Preconditioned_HUTCH(A,mi,eye(n,n),l,10);

        %stime finali
        HutchEst(l/10) = trPREC(l/10) + Hutch(l/10);
        HutchEstInv(l/10) = trPRECinv(l/10) + Hutchinv(l/10);

        PreHutchEst(l/10) = trPREC(l/10) + PreHutch(l/10);
        PreHutchEstInv(l/10) = trPRECinv(l/10) + PreHutchinv(l/10);

        BrutEst(l/10) = trPREC(l/10) + brutal_estimate(l/10);
        BrutEstInv(l/10) = trPRECinv(l/10) + brutal_estimate_inv(l/10);
        
    end
    
    figure(2)
    plot([10:10:10*dl],PreHutchEst,'r')
    hold on
    plot([10:10:10*dl],PreHutchEstInv,'-p')
    hold on
    plot([10:10:10*dl],HutchEst,'b')
    hold on
    plot([10:10:10*dl],HutchEstInv,'-v')
    hold on
    plot([10:10:10*dl],BrutEst,'g')
    hold on
    plot([10:10:10*dl],BrutEstInv,'y')
    hold on
    plot([10:10:10*dl],HutchA,'o-')
    hold on
    plot([10:10:10*dl],trA*ones(dl,1),'c')
    legend('PreHutch++','PreHutchInv++','Hutch++','HutchInv++','BrutalEst','BrutalEstInv','Hutch(A)','tr$(A+\mu I)$','interpreter','latex')
    
    figure(3)
    plot([10:10:10*dl],PreHutchEst/trA,'r')
    hold on
    plot([10:10:10*dl],PreHutchEstInv/trA,'-p')
    hold on
    plot([10:10:10*dl],HutchEst/trA,'b')
    hold on
    plot([10:10:10*dl],HutchEstInv/trA,'-v')
    hold on
    plot([10:10:10*dl],BrutEst/trA,'g')
    hold on
    plot([10:10:10*dl],BrutEstInv/trA,'y')
    hold on
    plot([10:10:10*dl],HutchA/trA,'o-')
    hold on
    ylim([0.8,1.2])
    legend('ratio PreHutch++','ratio PreHutchInv++','ratio Hutch++','ratio HutchInv++','ratio BrutalEst','ratio BrutalEstInv','ratio Hutch(A)')
end