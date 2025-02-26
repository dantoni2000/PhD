%test 21 02 2025
%on best error E: the parallel one! 

err=1e-03;
n=1000;
C = [ones(n,1); -ones(n,1)];
Ej = err*ones(2*n,1);
A = C*C';
%A = C*randn(2*n,1)';

CE=C+Ej;
E=zeros(2*n,2*n);
E(:,1)=Ej; E(:,2:2*n)=err*CE*ones(2*n-1,1)';

nc = norm(C)^2;
ne = norm(E,'fro')^2;

%controllo 1
ncce = norm(CE*pinv(CE)*E)^2;

%controllo 2
nf = norm(Ej*pinv(C)*A)^2;
npf = norm(Ej)^2*norm(pinv(C)*A)^2;

%controllo finale :(
n1 = norm(A-(CE)*pinv(CE)*(A+E),'fro')^2;
n2 = npf + ne;