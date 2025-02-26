%test 13 02 2025
%on column selection

n = 20;

u1 = ones(n,10);
u2 = 2*ones(n,9);
u3 = zeros(n,1);

A = [u1 u2 u3];

[U,S,V] = svd(A);

tA = A-A(:,10)*pinv(A(:,10))*A;
t2A = A-A(:,2)*pinv(A(:,2))*A;

[tU,tS,tV] = svd(tA);

