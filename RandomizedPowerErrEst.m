function E=RandomizedPowerErrEst(A,U,Lambda,q)
tol = 1e-3;
[n,~] = size(A);
w = randn(n,1);
y_0 = w/norm(w);
y = y_0;
csi(1) = 0;

for i = 1:q
    y_A = A*y; 
    y_U = U'*y;
    csi(i+1) = y'*y_A - y_U'*(diag(Lambda).*y_U);

    y = y_A - U*(diag(Lambda).*y_U);
    y = y/norm(y);

    if abs(csi(i+1) - csi(i)) < tol*csi(i)
       break
    end
end 

E=csi(end);