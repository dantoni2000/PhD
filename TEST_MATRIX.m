for r=1:10^4

c1 = 1e-1*randn; c2= 1e-1*randn; c=c1^2 + c2^2;
C = [c1; c2; 0];
err = 1e-4;
E = [-err^2/c * c1; -err^2/c * c2; sqrt(err^2 -err^4/c^2*c1^2 - -err^4/c^2*c2^2)];

Ct = C+E;
[Q,~] = qr(Ct,0);

n1 = norm(Q*Q'*C)^2;
n2 = norm(C)^2 - norm(E)^2;

dif(r) = n1 - n2;

end

max(dif)