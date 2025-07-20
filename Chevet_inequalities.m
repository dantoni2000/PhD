clear all
close all

for TT = 1:100

    n = 100; k = 20; T = 1000; q = 4;
    
    A = diag(randn(n,1)).^2;
    B = diag(randn(k,1)).^2;
    TrueNorm = 0;
    
    for t = 1:T
        Omega = randn(n,k);
        TrueNorm = TrueNorm + norm(A*Omega*B)^q;
    end
    
    TrueExp = 1/T*TrueNorm;
    Chev = (norm(A,'fro')*norm(B)+norm(A)*norm(B,'fro'))^q;
    
    r(TT) = Chev / TrueExp;

end

min(r)