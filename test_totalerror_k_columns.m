%test 21 02 2025
%on best error E: the parallel one! [Works well in practice] 

err=1e-04;
n=1000;

C1 = [ones(n,1); -ones(n,1)];
C2 = 2*C1;
C3 = 3*C1;
C = [C1 C2 C3];

E1 = err*ones(2*n,1);
E2 = 2*E1;
E3 = 3*E1;
EJ = [E1 E2 E3];

A = C*C';
%A = C*randn(2*n,1)';

CE=C+EJ;
E=err*CE*ones(2*n,3)'; E(:,1:3)=EJ; 

nc = norm(C,'fro')^2;
ne = norm(E,'fro')^2;

%controllo 1
ncce = norm(CE*pinv(CE)*E,'fro')^2;

%controllo 2
nf = norm(EJ*pinv(C)*A,'fro')^2;
npf = norm(EJ)^2*norm(pinv(C)*A,'fro')^2;

%controllo finale :(
n1 = norm(A-(CE)*pinv(CE)*(A+E),'fro')^2;
n2 = npf + ne;