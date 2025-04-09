function [its,Ptr] = Preconditioned_HUTCH(A,mi,P,N,m)

[n,~] = size(A);
Ptr = 0;
Pit = zeros(N,1);

for i = 1:N
    v = randsrc(n,1);
    % v = randn(n,1);
    [Pit(i),tr] = Preconditioned_Lanczos_log(A,mi,P,v,m);
    Ptr = Ptr + tr;
end

its = 1/N*sum(Pit);
Ptr=1/N*Ptr;

