clear all
close all

n = 300; T = 100;

for k = 20:20:280

    for t = 1:T
        Omega = randn(n,k);
        pinvOmega = pinv(Omega(1:k,1:end));
        TnpinvOmega(k/20,t) = norm(pinvOmega,'fro')^2;
        TSpecnpinvOmega(k/20,t) = norm(pinvOmega)^2;
        rap2(k/20,t) = TnpinvOmega(k/20,t)/k^2;
        rap1_5(k/20,t) = TnpinvOmega(k/20,t)/(5*k);
        rap1(k/20,t) = TnpinvOmega(k/20,t)/k;
        % FroSpec(k/900,t) = TnpinvOmega(k/900,t)/TSpecnpinvOmega(k/900,t);

    end
    
end

figure(1)
semilogy(rap2(12,:))
rap2(rap2(12,:)>1e+2)=1;
semilogy(rap2(12,:))
hold on
semilogy(1/T*sum(rap2(12,1:end))*ones(1,T))


figure(2)
semilogy(rap1_5(12,:))
rap2log(rap1_5(12,:)>1e+2)=1;
semilogy(rap1_5(12,:))
hold on
semilogy(1/T*sum(rap1_5(12,1:end))*ones(1,T))


figure(3)
semilogy(rap1(12,:))
rap1(rap1(12,:)>1e+2)=1;
semilogy(rap1(12,:))
hold on
semilogy(1/T*sum(rap1(12,1:end))*ones(1,T))


% figure(4)
% semilogy(FroSpec(1,:))
% hold on
% semilogy(1/T*sum(FroSpec(1,1:end))*ones(1,T))