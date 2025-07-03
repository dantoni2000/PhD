function [its,tr,h] = EST(A,mi,N,m,flag2)

if flag2==0
    [n,~] = size(A);
    tr = 0;
    it = zeros(N,1);
    
    for i = 1:N
        v = randsrc(n,1);
        % v = randn(n,1);
        [it(i),trr] = Lanczos_log(A,mi,v,m);
        tr = tr + trr(end,1);
        truetr = v'*logm(A+mi*eye(n))*v;
        abs_err = abs(trr-truetr);
        rel_err = abs_err/abs(truetr);
        h = semilogy(abs_err,'g');
        hold on
    end
    
    its = 1/N*sum(it);
    tr=1/N*tr;
    hold off

else
    [n,~] = size(A);
    tr = 0;
    it = zeros(N,1);
    
    for i = 1:N
        v = randsrc(n,1);
        %v = randn(n,1);
        [it(i),trr] = Lanczos_log(A,mi,v,m);
        tr = tr + trr(end,1);
    end
    
    its = 1/N*sum(it);
    tr=1/N*tr;
end