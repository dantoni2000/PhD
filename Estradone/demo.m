function demo()
    rng(1); % per riproducibilità
    n = 100;
    p = 0.05;
    A = gen_gnp(n,p);
    fprintf('Generato G(%d, %.3f). Numero spigoli = %d\n', n, p, nnz(triu(A,1)));
    fprintf('Estrada index (Gnp): %.6f\n', estrada_index(A));
    
    % Esempio multipartito
    parts = [25, 25, 50]; % somme = n
    p_between = 0.08;
    B = gen_random_multipartite(parts, p_between);
    fprintf('Generato grafico multipartito. Numero spigoli = %d\n', nnz(triu(B,1)));
    fprintf('Estrada index (multipartite): %.6f\n', estrada_index(B));
end