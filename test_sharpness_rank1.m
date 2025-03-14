%test 12 03 2025
%on the sharpness of the bound

clear all, close all

err = 1e-5;
n = 2;
s = 1;
delta=1e-4;
E_right=[];

T=1e+6;

factor=zeros(T,1);

V1 = [sqrt(delta+1/4*delta^2); 1-1/2*delta];
% V2 = [0; -sqrt(delta+1/4*delta^2); zeros(n-3,1); 1-1/2*delta];
Vs = V1';
Us=randn(n,1);
Us=Us/norm(Us);

A = randn(n,1)*V1';
na=norm(A,'fro')^2;

C = A(:,1:s);

VJ = pinv(V(1:s,1:s));

for rr = 1:T

    EJ = randn(n,s);
    EJ = sqrt(s)*err*EJ/norm(EJ,'fro');
    nej = s*err^2;

    [Q,~]=qr(C+EJ,0);

    
    %controllo 2
    nf = norm(A-Q*Q'*A,'fro')^2;
    npf = min(na, nej*norm(VJ,'fro')^2);
    factor(rr) = npf/nf;

    if factor(rr)<1.2
        E_right=[E_right EJ];
    end

end

best=min(factor);