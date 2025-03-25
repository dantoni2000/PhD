function tr = Lanczos_log(A,x,m)

u(:,1) = x/norm(x);
w(:,1) = A*u(:,1);
alpha(1) = u(:,1)'*w(:,1);
r(:,1) = w(:,1)-alpha(1)*u(:,1);
beta(1) = norm(r(:,1));
u(:,2) = r(:,1)/beta(1);


for i=2:m

    w(:,i) = A*u(:,i);
    alpha(i) = u(:,i)'*w(:,i);
    r(:,i) = w(:,i)-alpha(i)*u(:,i)-beta(i-1)*u(:,i-1);
    beta(i) = norm(r(:,i));
    u(:,i+1) = r(:,i)/beta(i);

end

T = spdiags([[beta(:,2:m); 0],alpha,beta],-1:1,m:m);
fT = log(T);
tr = norm(x)^2 * fT(1,1);