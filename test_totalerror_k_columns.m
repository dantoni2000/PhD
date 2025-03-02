%test 21 02 2025
%on best error E: the parallel one! [Works well in practice] 

clear all, close all
err=1e-04;
n=1000;

C1 = [ones(n,1); -ones(n,1)];
E1 = err*ones(2*n,1);

for j=1:3

C(:,j) = j*C1;
EJ(:,j) = j*E1;

end

A = C*C';
%A = C*randn(2*n,1)';

CE=C+EJ;
E=err*CE*ones(2*n,j)'; E(:,1:j)=EJ; 

nc = norm(C,'fro')^2;
ne = norm(E,'fro')^2;

%controllo 1
ncce = norm(CE*pinv(CE)*E,'fro')^2;

%controllo 2
nf = norm(EJ*pinv(C)*A,'fro')^2;
npf = norm(EJ,'fro')^2*norm(pinv(C)*A,'fro')^2;

%controllo finale :)
n1 = norm(A-(CE)*pinv(CE)*(A+E),'fro')^2;
n2 = npf + ne;