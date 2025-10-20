%% GRAFO MULTIPARTITO CASUALE
function A = gen_random_multipartite(parts, p_between, varargin)
% A = gen_random_multipartite(parts, p_between)
% parts: vettore dei tagli delle partizioni (es [n1 n2 ... nm])
% p_between: probabilità che esista un arco tra vertici in partizioni diverse
% Opzionale: p_within (default 0) probabilità di archi all'interno della stessa partizione
    ip = inputParser;
    addRequired(ip,'parts',@isvector);
    addRequired(ip,'p_between',@isnumeric);
    addParameter(ip,'p_within',0,@isnumeric);
    parse(ip,parts,p_between,varargin{:});
    p_within = ip.Results.p_within;

    parts = round(parts(:))';
    n = sum(parts);
    A = zeros(n);
    idx = [0 cumsum(parts)];
    for i = 1:length(parts)
        for j = i:length(parts)
            % indici nei blocchi
            rows = (idx(i)+1):idx(i+1);
            cols = (idx(j)+1):idx(j+1);
            if i == j
                % blocco diagonale: generiamo solo se p_within>0
                if p_within > 0
                    U = triu(rand(length(rows)),1);
                    T = double(U < p_within);
                    A(rows,cols) = T + T';
                end
            else
                % tra partizioni: probabilità p_between
                U = rand(length(rows), length(cols));
                T = double(U < p_between);
                A(rows,cols) = T;
                A(cols,rows) = T';
            end
        end
    end
end