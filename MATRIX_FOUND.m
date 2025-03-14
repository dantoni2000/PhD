% c1 = 1e-3; c2= 1e-2; c=c1^2 + c2^2;
% C = [c1; c2; 0];
% err = 1e-4;
% E = [-err^2/c * c1; -err^2/c * c2; sqrt(err^2 -err^4/c^2*c1^2 - -err^4/c^2*c2^2)];
% 
% Ct = C+E;
% [Q,~] = qr(Ct,0);
% 
% n1 = norm(Q*Q'*C)^2;
% n2 = norm(C)^2 - norm(E)^2;
% 
% n1 - n2

T=1e+4;

E0=[];

for s=1:T

    c1 = 1e-3; c2= 1e-3; c=c1^2 + c2^2;
    C = [c1; 3*c2];
    err = 1e-4;
    % E = [-err^2/c * c1; sqrt(err^2-err^4/c^2 * c1^2)];
    E = randn(2,1);
    E = E/norm(E)*err;
    
    Ct = C+E;
    [Q,~] = qr(Ct,0);
    
    n1 = norm(Q*Q'*C)^2;
    n2 = norm(C)^2 - norm(E)^2;
    
    m(s)=n1 - n2;
  
    if m(s)<1e-14
        E0=[E0 E];
    end

end