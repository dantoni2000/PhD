clear all
close all

n = 100; k = 20; maxp = 10; T = 1000;

expct = zeros(maxp,1);

for p = 1:maxp

    for t = 1:T
        Omega = randn(n,k+p-1);
        pinvOmega = pinv(Omega(1:k,1:end));
        SpecpinvOmega(p,t) = norm(pinvOmega)^4;
        TnpinvOmega(p,t) = norm(pinvOmega,'fro')^2;
        FourthTnpinvOmega(p,t) = norm(pinvOmega,'fro')^4;
        S4TnpinvOmega(p,t) = norm(pinvOmega'*pinvOmega,'fro')^2;
    end

    SpecnpinvOmega(p,1) = 1/T * sum(SpecpinvOmega(p,1:end));
    npinvOmega(p,1) = 1/T * sum(TnpinvOmega(p,1:end));
    FourthnpinvOmega(p,1) = 1/T * sum(FourthTnpinvOmega(p,1:end));
    S4npinvOmega(p,1) = 1/T * sum(S4TnpinvOmega(p,1:end));
    
    if p>2
        expct(p,1) = k/(p-2);
    end

    if p>4
        Fourthexpct(p,1) = sqrt(exp(1)^5/2)*(k/(p))^2;
        S4expct(p,1) = (k*(k+p-2))/((p-1)*(p-2)*(p-4));
    end
end