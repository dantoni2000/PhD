function [its,Ptr,h] = Preconditioned_HUTCH(A,mi,P,N,m,flag1,flag2)

if flag2==0
    [n,~] = size(A);
    Ptr = 0;
    Pit = zeros(N,1);
    
    for i = 1:N
        v = randsrc(n,1);
        % v = randn(n,1);
        [Pit(i),tr] = Preconditioned_Lanczos_log(A,mi,P,v,m);
        Ptr = Ptr + tr(end,1);
        truetr = v'*logm(P*(A+mi*eye(n))*P)*v;
        abs_err = abs(tr-truetr);
        rel_err = abs_err/abs(truetr);
    
        if flag1==0
            h = semilogy(abs_err,'r');
        elseif flag1==1
            h = semilogy(abs_err,'c');
        else
            h = semilogy(abs_err,'b');
        end
    
        hold on
    end
    
    its = 1/N*sum(Pit);
    Ptr=1/N*Ptr;
    hold off

else
    [n,~] = size(A);
    Ptr = 0;
    Pit = zeros(N,1);
    
    for i = 1:N
        v = randsrc(n,1);
        % v = randn(n,1);
        [Pit(i),tr] = Preconditioned_Lanczos_log(A,mi,P,v,m);
        Ptr = Ptr + tr(end,1);
        truetr = v'*logm(P*(A+mi*eye(n))*P)*v;
        abs_err = abs(tr-truetr);
        rel_err = abs_err/abs(truetr);
    end
    
    its = 1/N*sum(Pit);
    Ptr=1/N*Ptr;
end