function [its,tr,h] = EST(A,mi,N,m,flag2,tolrel,tolabs)

if nargin < 6 || isempty(tolrel)
        tolrel = 1e-10;   % valore di default per tolrel
end
if nargin < 7 || isempty(tolabs)
        tolabs = 1e-10;   % valore di default per tolabs
end

if flag2==0
    [n,~] = size(A);
    tr = 0;
    it = zeros(N,1);
    
    for i = 1:N
        % v = randsrc(n,1);
        v = randn(n,1);
        [it(i),trr] = Lanczos_log(A,mi,v,m,tolrel,tolabs);
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
        % v = randsrc(n,1);
        v = randn(n,1);
        [it(i),trr] = Lanczos_log(A,mi,v,m,tolrel,tolabs);
        tr = tr + trr(end,1);
    end
    
    its = 1/N*sum(it);
    tr=1/N*tr;
end