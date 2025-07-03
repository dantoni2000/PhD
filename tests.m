clear all
close all

n = 100; k = 10; maxp = 5; T = 100000;

expct = zeros(maxp,1);

for p = 1:maxp

    for t = 1:T
        Omega = randn(n,k+p-1);
        pinvOmega = pinv(Omega(1:k,1:end));
        TnpinvOmega(p,t) = norm(pinvOmega,'fro')^2;
    end

    npinvOmega(p,1) = 1/T * sum(TnpinvOmega(p,1:end));
    
    if p>2
        expct(p,1) = k/(p-2);
    end
end