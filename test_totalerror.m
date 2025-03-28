%test 21 02 2025
%on best error E: the parallel one! 

err=1e-02;
n=1000;
C = [ones(n,1); -ones(n,1)];
Ej = err*ones(2*n,1);
A = C*C';
%A = C*randn(2*n,1)';

CE=C+Ej;
E=zeros(2*n,2*n);
E(:,1)=Ej; E(:,2:2*n)=err*CE*ones(2*n-1,1)';
E=E*norm(Ej)/norm(E(:,2));

nc = norm(C)^2;
ne = norm(E,'fro')^2;

%controllo 1
ncce = norm(CE*pinv(CE)*E,'fro')^2;

%controllo 2
nf = norm(Ej*pinv(C)*A)^2;
npf = norm(Ej)^2*norm(pinv(C)*A,'fro')^2;

%controllo finale :)
n1 = norm(A-(CE)*pinv(CE)*(A+E),'fro')^2;

%n2 = npf*(nc^2 + norm(Ej)^2*nc)/(nc^2 + 2*norm(Ej)^2*nc +norm(Ej)^4) + (2*n-1)/(2*n)*ne + norm(Ej)^4*(nc + norm(Ej)^2)/(nc^2 + 2*norm(Ej)^2*nc +norm(Ej)^4);

n2= npf + ne;