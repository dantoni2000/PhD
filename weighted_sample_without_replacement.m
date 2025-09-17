function S = weighted_sample_without_replacement(prob, k)
    % prob: vettore di probabilità (non normalizzato)
    % k: numero di elementi da campionare
    n = length(prob);
    available = 1:n;
    probs_temp = prob(:)';  % vettore riga

    S = zeros(1, k);
    for t = 1:k
        cdf = cumsum(probs_temp);
        cdf = cdf / cdf(end);
        r = rand();
        idx = find(cdf >= r, 1, 'first');
        S(t) = available(idx);
        
        % rimuove l'elemento scelto
        available(idx) = [];
        probs_temp(idx) = [];
    end
end