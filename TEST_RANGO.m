%test 21 02 2025
%on best error E: for rank s it doesn't work well...

clear all, close all
err = 1e-02;
n = 1000;
s = 40;


CEJ = randn(1000,2*s);
[Q,R] = qr(CEJ,'econ');

% R = eye(s,s);
% R = triu(ones(s,s));
% R = triu(randn(s,s));

QC = Q(:,1:s);
QEJ = Q(:,s+1:2*s);

C = QC*R(1:s,1:s);
EJ = err*QEJ*R(1:s,1:s)/norm(R(1:s,1:s)'*R(1:s,1:s),'fro');

A = C*randn(n,s)'; A(:,1:s) = C;

CE = C+EJ;
E = err*CE*randn(n,s)'; E(:,1:s) = EJ;
E=E/norm(E,'fro')*norm(EJ,'fro');

nc = norm(C,'fro')^2;
ne = norm(E,'fro')^2;
nej = norm(EJ,'fro')^2;

%controllo 1
ncce = norm(CE*pinv(CE)*E,'fro')^2;

%controllo 2
nf = norm(EJ*pinv(C)*A,'fro')^2;
npf = norm(EJ,'fro')^2*norm(pinv(C)*A,'fro')^2;
npf/nf

%controllo finale :c
n1 = norm(A-(CE)*pinv(CE)*(A+E),'fro')^2;

n2 = npf*(nc^2 + nej*nc)/(nc^2 + 2*nej*nc +nej^2) + (2*n-s)/(2*n)*ne + nej^2*(nc + nej)/(nc^2 + 2*nej*nc +nej^2);


%n2 = npf + ne;

% RR=R(1:s,1:s)'*R(1:s,1:s);
% G=RR^-1 * trace(RR);
% eigs(G);