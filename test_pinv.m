%test 11 03 2025
%on best error E by choosing Ej as pseudoinverse of Vj...a bit better but
%not decisive yet!

clear all, close all
err = 1e-02;
n = 1000;
s = 100;

for rr = 1:100

CEJ = randn(n,2*s);

[Q,R] = qr(CEJ,'econ');

% R = eye(s,s);
% R = triu(ones(s,s));
% R = triu(randn(s,s));

QC = Q(:,1:s);
%QEJ = Q(:,s+1:2*s);

C = QC*R(1:s,1:s);

A = C*randn(n,s)'; A(:,1:s) = C;
[U,S,V]=svd(A);

S=S(1:s,1:s);
Vs=V(1:s,:);

VJ=pinv(C)*A*Vs';

EJ = U(:,s+1:2*s)*S*VJ';
EJ = sqrt(s)*err*EJ/norm(EJ,'fro');

CE = C+EJ;
E = CE*randn(n,s)'; E=sqrt(n)*err*E/norm(E,'fro'); E(:,1:s) = EJ;

nc = norm(C,'fro')^2;
ne = norm(E,'fro')^2;
nej = norm(EJ,'fro')^2;

%controllo 1
ncce = norm(CE*pinv(CE)*E,'fro')^2;

%controllo 2
nf = norm(EJ*pinv(C)*A,'fro')^2;
npf = norm(EJ,'fro')^2*norm(pinv(C)*A,'fro')^2;
factor(rr) = npf/nf;

%controllo finale :c
n1 = norm(A-(CE)*pinv(CE)*(A+E),'fro')^2;

n2 = npf*(nc^2 + nej*nc)/(nc^2 + 2*nej*nc +nej^2) + (2*n-s)/(2*n)*ne + nej^2*(nc + nej)/(nc^2 + 2*nej*nc +nej^2);

end

plot(1:rr, factor)
%n2 = npf + ne;

% RR=R(1:s,1:s)'*R(1:s,1:s);
% G=RR^-1 * trace(RR);
% eigs(G);