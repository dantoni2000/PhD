clear all
close all

n = 1000;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n)';
mi = 1;
G = 1000 * diag((-g+g(n,1)+1)./g(n,1));
% G = 1000 * diag(exp(-1*g));
A = Q*G*Q';

for k = 50:50:950
    [W,~] = qr(randn(n,k),'econ');
    Pj = W*W';
    
    M1=norm(logm(A+eye(n)) - sqrtm(logm(A+eye(n)))*Pj*sqrtm(logm(A+eye(n))),'fro');
    M2=norm(logm(A+eye(n)) - Pj*logm(A+eye(n))*Pj','fro');
    M3=norm(logm(A+eye(n)) - logm(Pj*A*Pj'+eye(n)),'fro');
    M4=norm(logm(A+eye(n)) - logm(sqrtm(A)*Pj*sqrtm(A)+eye(n)),'fro');

    r(k/50) = M4 / M1;
end 
plot(r)

sqrtm(logm(A+eye(n)))\logm(sqrtm(A)*Pj*sqrtm(A)+eye(n))/sqrtm(logm(A+eye(n))) - Pj