function [U,Lhat] = Adaptive_Nystrom(A,l_0,l_max,q,tol_e,tol_r,mi)

[n,~] = size(A);
Y = []; Om = []; E = 10^20; lambda_l = 10^20;
m = l_0;

while E > tol_e  &&  lambda_l/mi > tol_r

    Om_0 = randn(n,m);
    [Om_0,~] = qr(Om_0,0);
    Y_0 = A*Om_0;
    Y = [Y Y_0]; Om = [Om Om_0];
    nu = sqrt(n)*eps(norm(Y,2)); 
    Y_nu = Y + nu*Om;
    C = chol(Om'*Y_nu);
    B = Y_nu/C;
    [U,S,~] = svd(B,0);
    % Lhat = max(0,S^2-nu*eye(m,m));
    Lhat = max(0,S^2 - nu*eye(size(S)));
    lambda_l = Lhat(m,m);
    E = RandomizedPowerErrEst(A,U,Lhat,q);
    m = l_0; l_0 = 2*l_0;
    t = false;

    if l_0 > l_max
        l_0 = l_0 - m;
        m = l_max - l_0;
        % l_0 = l_max;
        % m = l_max;
        Om_0 = randn(n,m);
        [Om_0,~] = qr(Om_0,0);
        Y_0 = A*Om_0;
        Y = [Y Y_0]; Om = [Om Om_0];
        nu = sqrt(n)*eps(norm(Y,2)); 
        Y_nu = Y + nu*Om;
        C = chol(Om'*Y_nu);
        B = Y_nu/C;
        [U,S,~] = svd(B,0);
        % Lhat = max(0,S^2-nu*eye(m,m));
        Lhat = max(0,S^2 - nu*eye(size(S)));
        lambda_l = Lhat(m,m);
        t = true;
        break
    end

end
