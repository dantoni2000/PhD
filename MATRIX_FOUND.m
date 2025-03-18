%test 15 03 2025 
%for a matrix C, we show that there exists Ej so that the bound is tight

T = 1e+4;

E0 = [];

for s=1:T

    c1 = 1e-3; c2= 1e-3; c=c1^2 + c2^2;
    C = [c1; -0.5*c2; c1; c1^2];
    err = 1e-4;
    E = randn(4,1);
    E = E/norm(E)*err;
    
    Ct = C+E;
    [Q,~] = qr(Ct,0);
    
    n1 = norm(Q*Q'*C)^2;
    n2 = norm(C)^2 - norm(E)^2;
    
    m(s) =abs(n1/n2);
  
    if m(s)<1+1e-09
        E0 = [E0 E];
    end

end