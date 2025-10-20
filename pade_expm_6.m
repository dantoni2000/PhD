function X = pade_expm_6(A)
% PADE_EXPM_6  Matrix exponential using Pade [6/6] approximation
% and scaling & squaring
%
% Reference: Higham (2005), Alg. 5.2 simplified for m = 6

n = size(A,1);
I = eye(n);

% --- Coefficienti della serie di e^x fino al grado 6 ---
c = [1.0, ...
     0.5, ...
     1/12, ...
     1/60, ...
     1/360, ...
     1/2520, ...
     1/20160];

% --- Step 1: scaling ---
normA = norm(A,1);
theta6 = 3.0;    % soglia empirica per m=6 (Higham tabella 3.1)
s = max(0, ceil(log2(normA/theta6)));

B = A / 2^s;

% --- Step 2: potenze ---
B2 = B*B;
B4 = B2*B2;
B6 = B2*B4;

% --- Step 3: U e V ---
U = B*(c(2)*I + c(4)*B2 + c(6)*B4);
U = U + c(1)*B;  % parte lineare

V = c(1)*I + c(3)*B2 + c(5)*B4 + c(7)*B6;

% --- Step 4: risolvi (V - U) X = (V + U)
X = (V - U) \ (V + U);

% --- Step 5: squaring
for k = 1:s
    X = X * X;
end

end