
n = 1000;
mi = 1;
alpha = 1; nu = 7/2; % nu = 13/2;
kernel = @(x,y) sqrt(pi)*((alpha*norm(x-y))^(nu)*besselk(nu,alpha*norm(x-y)))/(2^(nu-1)*alpha^(2*nu)*gamma(nu+0.5));
data_matrix = randn(1,n);

for row = 1:n

    for column = 1:n

        if row == column
        
            A(row,column) = sqrt(pi)*gamma(nu)/(gamma(nu + 1/2)*alpha^(2*nu));

        else

            A(row,column) = kernel(data_matrix(:,row),data_matrix(:,column));

        end

    end

end

[Q,G] = svd(A);
G=10^10*diag(diag(G));
% G = G + .1*eye(n);
A = Q*G*Q'; 
% A = G;

n = 1000;
Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n)';
mi = 1;
% G = 10^8*diag(exp(-.12*g));
% G = 1000*diag(exp(-.1*g));
A = Q*G*Q';
% A = G;

% G = logm(G+eye(n));

j = 64;
% nuc = 0*sum(diag(G(j+6:end,j+6:end).^2)) + 1.5*sum(diag(G(j+7:end,j+7:end).^2));
% nuc2 = sum(diag(G(j+6:end,j+6:end).^2)) + 0.5*sum(diag(G(j+7:end,j+7:end).^2));

nuc3=sum(diag(G(j+6:end,j+6:end).^2));

for k = 1:n-j-6
    nuc3 = nuc3 + 2*sum(diag(G(j+6+k:end,j+6+k:end).^2));
end

nuc3 = nuc3/2;

fro = sum(diag(G(j:j+5,j:j+5).^2));


% fro - nuc
% fro / nuc
% 
% fro - nuc2
% fro/nuc2

fro - nuc3
fro/nuc3