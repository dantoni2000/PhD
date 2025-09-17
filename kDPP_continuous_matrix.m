%% 2. Funzione kDPP_continuous_matrix (continua)
function X = kDPP_continuous_matrix(K, k)
    n = size(K, 1);
    [V, D] = eig(K);
    lambda = diag(D);
    
    % Seleziona k autovettori pesati
    J = weighted_sample_without_replacement(lambda / sum(lambda), k);
    V_selected = V(:, J);
    
    % Campionamento iterativo subset
    S = zeros(1, k);
    Vt = V_selected;
    for t = 1:k
        probs = sum(Vt.^2, 2);
        probs = probs / sum(probs);
        cdf = cumsum(probs);
        cdf = cdf / cdf(end);
        r = rand();
        i = find(cdf >= r, 1, 'first');
        S(t) = i;
        
        % Aggiornamento Gram-Schmidt
        e = Vt(i, :)';
        Vt = Vt - (Vt * e) * e';
        Vt = orth(Vt);
    end
    
    % Costruzione matrice continua
    X = zeros(n, n);
    X(S, S) = K(S, S);
end