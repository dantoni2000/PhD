function [its,tr]=Preconditioned_Lanczos_log(A,mi,P,x,m)

tol=1e-12;
its = m;
tr(1,1) = 0; 

u(:,1) = x/norm(x);
Pu(:,1) = P*u(:,1);
w(:,1) = P*(A*Pu(:,1) + mi*Pu(:,1));
alpha(1) = u(:,1)'*w(:,1);
r(:,1) = w(:,1)-alpha(1)*u(:,1);
beta(1) = norm(r(:,1));
u(:,2) = r(:,1)/beta(1);
T(1,1) = alpha(1);

for i=2:m
    oldtr = tr(i-1,1);
    Pu(:,i) = P*u(:,i);
    w(:,i) = P*(A*Pu(:,i) + mi*Pu(:,i));
    alpha(i,1) = u(:,i)'*w(:,i);
    r(:,i) = w(:,i)-alpha(i)*u(:,i)-beta(i-1)*u(:,i-1);

    % Ortogonalizzazione con GS
    for j=1:i
        r(:,i) = r(:,i)-(u(:,j)'*r(:,i))*u(:,j);
    end

    % Riortogonalizzazione con GS
    for j=1:i
        r(:,i) = r(:,i)-(u(:,j)'*r(:,i))*u(:,j);
    end

    beta(i,1) = norm(r(:,i));
    u(:,i+1) = r(:,i)/beta(i);

    T(i,i) = alpha(i,1); T(i,i-1)=beta(i-1,1); T(i-1,i)=beta(i-1,1);
    fT = logm(T);
    tr(i,1) = norm(x)^2 * fT(1,1);

    if abs(tr(i,1)-oldtr)<=tol*abs(tr(i,1))
        its=i-1;
        break
    end

end
% T = diag(alpha) + diag(beta(1:end-1),1) + diag(beta(1:end-1),-1);
% fT = logm(T);
% tr = norm(x)^2 * fT(1,1);