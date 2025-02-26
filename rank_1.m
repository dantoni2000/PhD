%test 16 02 2025
%on rank 1 matrices with incoherent singular vectors

n=100;

case_ze = input('If case=0 zero matrices, If case=1 eps matrices: ');

switch case_ze
    case 0
U=zeros(n,1);
U(randperm(n,10))=1/sqrt(10);
S=20;
V=zeros(n,1);
V(randperm(n,10))=1/sqrt(10);
end 

switch case_ze
    case 1
rf=1e-10;
U=rf*ones(n,1);
U(randperm(n,10))=1/sqrt(10);
S=20;
V=rf*ones(n,1);
V(randperm(n,10))=1/sqrt(10);
end 

A=U*S*V';
pb=0;

for j=1:1000
    V_c=zeros(n,1);
    r=randperm(n,1);

    V_c(r)=V(r);
    A_c=U*S*V_c';

    if norm(A-A_c*pinv(A_c)*A)<1e-12
        pb=pb+1;
    end
end

