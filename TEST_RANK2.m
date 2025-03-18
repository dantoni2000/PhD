T=1e+6;

E0=[];

for s=1:T

    c1 = 1e-2; c2= 1e-2;
    C1 = [c1; 3*c2; c1];
    C2 = [-2*c2; c1; c1];
    C = [C1 C2];
    nc = norm(C,'fro')^2;
    err = 5e-3;
    % E1 = [-err^2/nc*c1; -err^2/nc*3*c2; sqrt(err^2*norm(C1)^2-(err^2/nc*c1)^2 -(err^2/nc*3*c2)^2); 0];
    % E2 = [+err^2/nc*2*c2; -err^2/nc*c1; 0; sqrt(err^2*norm(C2)^2-(err^2/nc*2*c2)^2 -(err^2/nc*c1)^2)];
    % E = [E1 E2];
    E = randn(3,2);
    E = E/norm(E)*err;
    
    Ct = C+E;
    [Q,~] = qr(Ct,0);
    
    n1 = norm(Q*Q'*C,'fro')^2;
    n2 = norm(C,'fro')^2 - norm(E,'fro')^2;
    
    m(s) = abs(n1/n2);
  
end

M=min(m)