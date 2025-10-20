%% MATRICE ADIACENZA MEDIA / CENTRALIZZATA (opzionale)
function B = centered_adjacency(A, p)
% B = centered_adjacency(A,p)
% Restituisce A - p*J (matrice con media p) utile per confronti con RMT.
% p è la probabilità di edge (se nota). J è matrice con zeri sulla diagonale.
    if nargin < 2
        error('Serve p (probabilità media) per centralizzare.');
    end
    n = size(A,1);
    J = ones(n) - eye(n);
    B = A - p*J;
end