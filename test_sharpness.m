%test 12 03 2025
%on the sharpness of the bound

clear all, close all

err = 1e-4;
n = 4;
s = 2;
delta=1e-6;
E_right=[];

T=1e+6;

factor=zeros(T,1);

[U,~]=qr(randn(n,s),0);

V1 = [sqrt(delta+1/4*delta^2); 0; 1-1/2*delta; 0];
% V1 = [1-1/2*delta; 0; sqrt(delta+1/4*delta^2); zeros(n-3,1)];
V2 = [0; sqrt(delta+1/4*delta^2); 0; 1-1/2*delta];
% V2 = [0;  1-1/2*delta; zeros(n-3,1);-sqrt(delta+1/4*delta^2)];

Vs = [V1 V2]';

A = U*[V1 V2]';

C = A(:,1:s);
na = norm(A,'fro')^2;

VJ = pinv(C)*A*Vs';

for rr = 1:T

EJ = rand(n,s);
EJ = sqrt(s)*err*EJ/norm(EJ,'fro');
nej = norm(EJ)^2;

%controllo 
nf = norm((eye(n,n)-(C+EJ)*pinv(C+EJ))*C*pinv(C)*A,'fro')^2;
npf = nej*norm(pinv(C)*A,'fro')^2;
factor(rr) = npf/nf;

% if factor(rr)<1.2
%     E_right=[E_right EJ];
% end

end

best=min(factor);