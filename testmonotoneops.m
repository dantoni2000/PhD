rng(42); % Fissa il seed per riproducibilità
n_trials = 1000000;
best_ratio = -Inf;

for i = 1:n_trials
    % Genera matrici ortogonali casuali
    [Q1, ~] = qr(randn(3));
    [Q2, ~] = qr(randn(3));
    
    % Autovalori positivi casuali in [1, 3]
    D1 = diag(1 + 2 * rand(3, 1));
    D2 = diag(1 + 2 * rand(3, 1));
    
    % Costruisci A, B SPD
    A = Q1 * D1 * Q1';
    B = Q2 * D2 * Q2';
    D = A - B;
    
    % Verifica che A - B sia SPD
    if issymmetric(D) && all(eig(D) > 0)
        fA = logm(eye(3) + A);
        fB = logm(eye(3) + B);
        
        lhs = norm(fA - fB, 'fro');
        rhs = norm(A - B, 'fro');
        ratio = lhs / rhs;
        
        if ratio > best_ratio
            best_ratio = ratio;
            best_A = A;
            best_B = B;
            best_fA = fA;
            best_fB = fB;
            best_lhs = lhs;
            best_rhs = rhs;
        end
    end
end

% Stampa il risultato più vicino alla violazione
fprintf('Miglior rapporto f(A)-f(B)/A-B: %.6f\n', best_ratio);
fprintf('lhs = %.6f\n', best_lhs);
fprintf('rhs = %.6f\n', best_rhs);
disp('Matrice A:'); disp(best_A);
disp('Matrice B:'); disp(best_B);