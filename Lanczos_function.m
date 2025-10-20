function [its,tr] = Lanczos_function(A,mi,x,m,tolrel,tolabs,fhandle)
its = m;
tr(1,1) = 0; 

u(:,1) = x/norm(x);
w(:,1) = A*u(:,1) + mi*u(:,1);
alpha(1) = u(:,1)'*w(:,1);
r(:,1) = w(:,1)-alpha(1)*u(:,1);
beta(1) = norm(r(:,1));
u(:,2) = r(:,1)/beta(1);
T(1,1) = alpha(1);

for i=2:m
    oldtr = tr(i-1,1);
    w(:,i) = A*u(:,i) + mi*u(:,i);
    alpha(i,1) = u(:,i)'*w(:,i);
    r(:,i) = w(:,i)-alpha(i)*u(:,i)-beta(i-1)*u(:,i-1);

    % Riortogonalizzazione con GS
    for j=1:i
        r(:,i) = r(:,i)-(u(:,j)'* r(:,i))*u(:,j);
    end

    beta(i,1) = norm(r(:,i));
    u(:,i+1) = r(:,i)/beta(i);

    T(i,i) = alpha(i,1); T(i,i-1)=beta(i-1,1); T(i-1,i)=beta(i-1,1);
    fT = fhandle(T);
    tr(i,1) = norm(x)^2 * fT(1,1);

    if abs(tr(i,1)-oldtr)<=tolrel*abs(tr(i,1))
        its=i-1;
        break
    end

    if abs(tr(i,1)-oldtr)<=tolabs
        its=i-1;
        break
    end

end

% T = diag(alpha) + diag(beta(1:end-1),1) + diag(beta(1:end-1),-1);
% fT = logm(T);
% tr = norm(x)^2 * fT(1,1);