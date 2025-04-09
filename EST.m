function [its,tr] = EST(A,mi,N,m)

[n,~] = size(A);
tr = 0;
it = zeros(N,1);

for i = 1:N
    v = randsrc(n,1);
    % v = randn(n,1);
    [it(i),trr] = Lanczos_log(A,mi,v,m);
    tr = tr + trr;
end

its = 1/N*sum(it);
tr=1/N*tr;