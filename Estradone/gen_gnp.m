%% ERDOS-RENYI G(n,p)
function A = gen_gnp(n, p, varargin)
% A = gen_gnp(n,p)
% Genera una matrice di adiacenza nxn (simmetrica, 0 diag) per modello G(n,p)
% Opzionale: 'directed',true/false (default false)
    ip = inputParser;
    addRequired(ip,'n',@isnumeric);
    addRequired(ip,'p',@isnumeric);
    addParameter(ip,'directed',false,@islogical);
    parse(ip,n,p,varargin{:});
    directed = ip.Results.directed;

    % matrice superiore casuale
    U = rand(n);
    if directed
        A = double(U < p);    % include possibili auto-archi, rimuoveremo diag sotto
        A(1:n+1:end) = 0;     % no loops
    else
        U = triu(U,1);
        T = double(U < p);
        A = T + T';           % simmetrica
    end
end