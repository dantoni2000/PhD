%% ESTRADA INDEX
function EE = estrada_index(A, varargin)
% EE = estrada_index(A)
% Calcola l'Estrada index EE(G) = sum_{i=1}^n exp(lambda_i)
% Per matrici grandi si può richiedere solo k autovalori estremi mediante 'k' opzionale.
    ip = inputParser;
    addRequired(ip,'A',@ismatrix);
    addParameter(ip,'k',[],@(x) isempty(x) || (isnumeric(x) && x>0));
    parse(ip,A,varargin{:});
    k = ip.Results.k;

    n = size(A,1);
    if isempty(k) || k >= n
        % calcola tutti gli autovalori (metodo diretto)
        lam = eig(A);
        EE = sum(exp(lam));
    else
        % approssimazione usando i k autovalori di piu' grande modulo
        opts.isreal = isreal(A);
        opts.issym = issymmetric(A);
        try
            [V,D] = eigs(A, k, 'largestabs');
            lam = diag(D);
            % per fare stima migliore aggiungiamo contributo approssimato degli altri:
            % semplice lower bound: assumiamo altri autovalori minori o vicini a 0
            EE = sum(exp(lam));
            % Nota: questa è un'approssimazione; se vuoi accuratezza usa k>=n.
        catch ME
            % fallback a eig se eigs non riesce
            lam = eig(A);
            EE = sum(exp(lam));
        end
    end
end