%test 21 02 2025
%on best error E: the parallel one! 

clear all, close all

err = 1e-02;
n = 1000;
s = 1;
k = 5;
CEJ = randn(1000,2*s);
[Q,~] = qr(CEJ,'econ');

% R = eye(s,s);
% R = triu(ones(s,s));
R = triu(randn(2*s,k));

QC = Q(:,1:s);
QEJ = Q(:,s+1:2*s);

C = QC*R(1:s,1:k);
EJ = sqrt(k)*err*QEJ*R(1:s,1:k)/norm(R(1:s,1:k),'fro');

%A = C*C';
A = C*randn(2*n,k)';

CE = C+EJ;
E = CE*ones(2*n,k)'; E(:,1:k)=EJ;
E = sqrt(n-k)*err*E/norm(E,'fro');

nc = norm(C,'fro')^2;
ne = norm(E,'fro')^2;
nej = norm(EJ,'fro')^2;


%controllo 1
ncce = norm(CE*pinv(CE)*E,'fro')^2;

%controllo 2
nf = norm(EJ*pinv(C)*A,'fro')^2;
npf = norm(EJ,'fro')^2*norm(pinv(C)*A,'fro')^2;

%controllo finale :)
n1 = norm(A-(CE)*pinv(CE)*(A+E),'fro')^2;

%n2 = npf*(nc^2 + nej*nc)/(nc^2 + 2*nej*nc +nej^2) + (2*n-3)/(2*n)*ne + nej^2*(nc + nej)/(nc^2 + 2*nej*nc +nej^2);

n2 = npf + ne;

%n1 - n2