%test 16 02 2025
%on rank (1 + noise) matrices with incoherent singular vectors

n=100;
ee=1e-10;
U=ee*ones(n,1);
U(randperm(n,10))=1/sqrt(10);
S=20;
V=ee*ones(n,1);
V(randperm(n,10))=1/sqrt(10);

rf=1e-14;
E=randn(n,n);
E=E/norm(E);
E=rf*E;


A=U*S*V';
AE=A+E;
[U,S,V]=svd(AE);
pb=0;

for j=1:1000
    V_c=zeros(n,n);
    r=randperm(n,1);

    V_c(r,:)=V(r,:);
    AE_c=U*S*V_c';

    if norm(AE-AE_c*pinv(AE_c)*AE)<1e-12
        pb=pb+1;
    end
end

