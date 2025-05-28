function [its,tr]=Nystrom_Lanczos_log(A,mi,U,L,x,m)

tol_beta = 1e-14;
tol=1e-12;
its = m;
tr(1,1) = 0; 

L = diag(diag(L+mi*eye(size(L))).^-0.5);

u(:,1) = x/norm(x);
Pu(:,1) = U*(diag(L).*(U'*u(:,1))) + mi * (u(:,1) - U*(U'*u(:,1)));
APu(:,1) = A*Pu(:,1) + mi*Pu(:,1);
w(:,1) = U*(diag(L).*(U'*APu(:,1))) + mi * (APu(:,1) - U*(U'*APu(:,1)));
alpha(1) = u(:,1)'*w(:,1);
r(:,1) = w(:,1)-alpha(1)*u(:,1);
beta(1) = norm(r(:,1));
if beta(1) < tol_beta
    tr(1,1) = norm(x)^2 * log(alpha(1));
    its = 1;
    return
end
u(:,2) = r(:,1)/beta(1);
T(1,1) = alpha(1);

for i=2:m
    oldtr = tr(i-1,1);
    Pu(:,i) = U*(diag(L).*(U'*u(:,i))) + mi * (u(:,i) - U*(U'*u(:,i)));
    APu(:,i) = A*Pu(:,i) + mi*Pu(:,i);
    w(:,i) = U*(diag(L).*(U'*APu(:,i))) + mi * (APu(:,i) - U*(U'*APu(:,i)));
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
      if beta(i) < tol_beta
        T(i,i) = alpha(i); 
        T(i,i-1) = beta(i-1); 
        T(i-1,i) = beta(i-1);
        fT = logm(T(1:i, 1:i));
        tr(i,1) = norm(x)^2 * fT(1,1);
        its = i;
      break
    end
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