function [its,Ptr,h] = Nystrom_HUTCH(A,mi,s,N,m,flag1,flag2)

if flag2==0
    [n,~] = size(A);
    Ptr = 0;
    Pit = zeros(N,1);

    [U,Lhat] = PinvNystrom(A,s);
    P = (mi)^0.5*U*(Lhat+mi*eye(s))^-0.5*U' + (eye(n) - U*U');

    Pp = sum(log(diag(Lhat+mi*eye(s,s))),"all");
    % lgPAP = logm(P*(A+mi*eye(n))*P);
    
    for i = 1:N
        v = randsrc(n,1);
        % v = randn(n,1);
        [Pit(i),tr] = Nystrom_Lanczos_log(A,mi,U,Lhat,v,m);
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
    Ptr = 1/N*Ptr + Pp;
    hold off

else
    [n,~] = size(A);
    Ptr = 0;
    Pit = zeros(N,1);

    [U,Lhat] = PinvNystrom(A,s);
    Pp = sum(log(diag(Lhat+mi*eye(s,s))),"all");
    
    for i = 1:N
        v = randsrc(n,1);
        % v = randn(n,1);
        [Pit(i),tr] = Nystrom_Lanczos_log(A,mi,U,Lhat,v,m);
        % P = (mi)^0.5*U*(Lhat+mi*eye(s))^-0.5*U' + (eye(n) - U*U');
        % tr(end,1) - trace(logm(P*(A+eye(n))*P))
        % Pp - trace(logm(A+eye(n))) + trace(logm(P*(A+eye(n))*P))
        % -2*trace(logm(P)) - trace(logm(A+eye(n))) + trace(logm(P*(A+eye(n))*P))
        Ptr = Ptr + tr(end,1);
    end
    
    its = 1/N * sum(Pit);
    Ptr = 1/N * Ptr + Pp;

end