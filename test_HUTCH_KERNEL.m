% 30 03 2025 Test Hutchinson: in quali ipotesi lo stimatore gratuito funziona meglio
% DATI KERNEL

clear all
close all

n = 1000;
dl = 10;
mi = 1e-2;

kernel = @(x,y) ((norm(x)/n)^10 + (norm(y)/n)^10)^(0.1);
data_matrix = orth(randn(n,n))*diag(1:n);
A = build_kernel_matrix(data_matrix,kernel);
[Q,G] = svd(A);

A = Q*G*Q';
figure(1)
semilogy(diag(G))
hold on 
semilogy(mi*ones(n))

trA = sum(log(diag(G+mi*eye(n,n))),"all");

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
    brutal_estimate(l/10) = n*log(mi);
    
    PreHutchEst(l/10) = trPREC(l/10) + PreHutch(l/10);
    HutchEst(l/10) = trPREC(l/10) + Hutch(l/10);
    BrutEst(l/10) = trPREC(l/10) + brutal_estimate(l/10);
    
end

figure(2)
plot([10:10:10*dl],PreHutchEst)
hold on
plot([10:10:10*dl],HutchEst)
hold on
plot([10:10:10*dl],BrutEst)
hold on
plot([10:10:10*dl],trA*ones(dl,1))
legend('PreHutch++','Hutch++','BrutalEst','tr$(A+\mu I)$','interpreter','latex')

figure(3)
plot([10:10:10*dl],PreHutchEst/trA)
hold on
plot([10:10:10*dl],HutchEst/trA)
hold on
plot([10:10:10*dl],BrutEst/trA)
ylim([0.8,1.2])
legend('ratio PreHutch++','ratio Hutch++','ratio BrutalEst')