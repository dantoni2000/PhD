function tr = Lanczos_log(A,x,m)

u(:,1) = x/norm(x);
w(:,1) = A*u(:,1);
alpha(1) = u(:,1)'*w(:,1);
r(:,1) = w(:,1)-alpha(1)*u(:,1);
beta(1) = norm(r(:,1));
u(:,2) = r(:,1)/beta(1);


for i=2:m

    w(:,i) = A*u(:,i);
    alpha(i,1) = u(:,i)'*w(:,i);
    r(:,i) = w(:,i)-alpha(i)*u(:,i)-beta(i-1)*u(:,i-1);
    beta(i,1) = norm(r(:,i));
    u(:,i+1) = r(:,i)/beta(i);

end
T = diag(alpha) + diag(beta(1:m-1),1) + diag(beta(1:m-1),-1);
fT = logm(T);
tr = norm(x)^2 * fT(1,1);