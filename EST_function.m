function [its,tr,h] = EST_function(A,mi,N,m,flag2,fhandle,tolrel,tolabs)

if nargin < 7 || isempty(tolrel)
        tolrel = 1e-10;   % valore di default per tolrel
end
if nargin < 8 || isempty(tolabs)
        tolabs = 1e-10;   % valore di default per tolabs
end

if flag2==0
    [n,~] = size(A);
    tr = 0;
    it = zeros(N,1);
    
    for i = 1:N
        % v = randsrc(n,1);
        v = randn(n,1);
        [it(i),trr] = Lanczos_function(A,mi,v,m,tolrel,tolabs,fhandle);
        tr = tr + trr(end,1);
        truetr = v'*fhandle(A+mi*eye(n))*v;
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
        [it(i),trr] = Lanczos_function(A,mi,v,m,tolrel,tolabs,fhandle);
        tr = tr + trr(end,1);
    end
    
    its = 1/N*sum(it);
    tr=1/N*tr;
end