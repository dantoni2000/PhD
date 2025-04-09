% 30 03 2025 Test Hutchinson: in quali ipotesi lo stimatore gratuito funziona meglio
% DATI SINTETICI

clear all
close all

T = 10;
n = 500;
dl = 12;

Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n);
mi = 1e-4;

decadimento=1;

switch decadimento
    case 1
        G = diag(1./(g).^2);
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
        
    end

switch decadimento
    case 2
        G = diag(1./sqrt(g));
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
    end

switch decadimento
    case 3
        G = diag(exp(-g));
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
    end

A = Q*G*Q';
trA = sum(log(diag(G+mi*eye(n,n))),"all");

for t=1:T
    for l=10:10:10*dl
    
        [U,Lhat] = Nystrom(A,l); ll=Lhat(l,l);
        P05 = sqrtm((ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U') );
        P = (ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U');
    
        trPREC(l/10) = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l))));
    
        true_value(l/10) = trace(logm(P05*(A+mi*eye(n))*P05));
    
        Hutch(l/10) = 0;
        
        for j=1:l
            
            x=randn(n,1);
            Hutch(l/10) = Hutch(l/10) + Lanczos_log((P05*(A+mi*eye(n))*P05),x,3);
            
        end
        
        PreHutch(l/10) = Preconditioned_HUTCH(A,mi,P,l,3);
        Hutch(l/10)= 1/l*Hutch(l/10);
        brutal_estimate(l/10) = n*log(mi); %+G(n,n)
        
        PreHutchEst(l/10) = trPREC(l/10) + PreHutch(l/10);
        HutchEst(l/10) = trPREC(l/10) + Hutch(l/10);
        BrutEst(l/10) = trPREC(l/10) + brutal_estimate(l/10);
        
    end
    
    figure(2)
    plot([10:10:10*dl],PreHutchEst,'r')
    hold on
    plot([10:10:10*dl],HutchEst,'b')
    hold on
    plot([10:10:10*dl],BrutEst,'y')
    hold on
    plot([10:10:10*dl],trA*ones(dl,1),'c')
    legend('PreHutch++','Hutch++','BrutalEst','tr$(A+\mu I)$','interpreter','latex')
    
    figure(3)
    plot([10:10:10*dl],PreHutchEst/trA,'r')
    hold on
    plot([10:10:10*dl],HutchEst/trA,'b')
    hold on
    plot([10:10:10*dl],BrutEst/trA,'y')
    ylim([0.8,1.2])
    legend('ratio PreHutch++','ratio Hutch++','ratio BrutalEst')
end