delta = 0.01;
A = [1, 1, 1; 0, delta, 1; 0, 0, 0];
C = A(:, 1:2);
[Q, ~] = qr(C, 0);
CCdaggerA = Q*Q'*A;
[U, S, V] = svd(A, 0);
V = V(:,1:2);

eps = 1e-4;

normV1Jplus = norm(inv(V(1:2,:)), 'fro')^2;

maxratio = 0;
for t = 1:1000000
    E = randn(3, 2);
    E = E / norm(E, 'fro') * eps;
    Ctilde = C + E;
    [Qtilde, ~] = qr(Ctilde, 0);
    ratio = norm(CCdaggerA - Qtilde * Qtilde' * CCdaggerA, 'fro')^2 ...
        / normV1Jplus / eps^2;
    if ratio > maxratio
        disp(t)
        disp(ratio)
        disp(E)
        disp(Qtilde)
    end
    maxratio = max(maxratio, ratio);
end

disp(maxratio)