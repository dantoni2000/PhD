%% ---------- Funzione: costruisce matrice (k+1)x(k+1) ----------
function A = construct_simplex_matrix(k, eps_shift, colnorm)
    % A = construct_simplex_matrix(k, eps_shift, colnorm)
    % k           : dimensione del simplex (numero di colonne = k+1)
    % eps_shift   : piccolo valore per lo shift ortogonale (es. 1e-6 .. 1e-2)
    % colnorm     : norma desiderata per ciascuna colonna del simplex (default 1)
    %
    % Restituisce A di dimensione (k+1) x (k+1).
    
    if nargin < 3
        colnorm = 1;
    end
    n = k + 1;
    
    % Matrice V: colonne = vertici di un simplex centrato (sommano a zero)
    V = eye(n) - (1/n) * ones(n);
    % Normalizza colonne in modo che ciascuna abbia norma colnorm
    for j = 1:n
        v = V(:, j);
        v = v / norm(v) * colnorm;
        V(:, j) = v;
    end
    
    % Direzione ortogonale (la direzione dei vettori 1): unit vector
    w = ones(n, 1) / sqrt(n);
    
    % Shift ortogonale: aggiunge un piccolo componente lungo w a tutte le colonne
    A = V + eps_shift * (w * ones(1, n));
end