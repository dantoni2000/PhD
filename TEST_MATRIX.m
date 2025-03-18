%test 15 03 2025 
%for a matrix C having a zero element, we explicitly found Ej so that the bound is tight

T = 1e+6;
CC = [];

for r=1:T

    c1 = 1e-1*randn; c2= 1e-1*randn; c=c1^2 + c2^2;
    C = [c1; c2; 0];
    err = 1e-4;
    E = [-err^2/c * c1; -err^2/c * c2; sqrt(err^2 -err^4/c^2*c1^2 - -err^4/c^2*c2^2)];
    
    Ct = C+E;
    [Q,~] = qr(Ct,0);
    
    n1 = norm(Q*Q'*C)^2;
    n2 = norm(C)^2 - norm(E)^2;
    
    dif(r) = n1 - n2;
    
    if dif(r)>1e-9
        CC = [CC C];
    end

end

max(dif)