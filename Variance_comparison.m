clear all
close all
warning off

T = 1; dl = 20; s = 10;

decadimento=0;

switch decadimento
    case 0
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-2;
        G = diag((-g+g(n,1)+1)./g(n,1));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 1
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-3;
        G = diag(1./(g).^2);
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 2
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-2;
        G = diag(1./sqrt(g));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 3
        n = 1000;
        Q = randn(n,n); [Q,~] = qr(Q);
        g = linspace(1,n,n)';
        mi = 1e-3;
        G = diag(exp(-g));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 4
        n = 1000;
        mi = 1e-4;
        kernel = @(x,y) ((norm(x)/n)^10 + (norm(y)/n)^10)^(0.1);
        data_matrix = orth(randn(n,n))*diag(1:n);
        A = build_kernel_matrix(data_matrix,kernel);
        [Q,G] = svd(A);
        G=diag(diag(G));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))

    case 5
        n = 1000;
        mi = 5e-6;
        alpha = 1; nu = 5/2;
        kernel = @(x,y) sqrt(pi)*((alpha*norm(x-y))^(nu)*besselk(nu,alpha*norm(x-y)))/(2^(nu-1)*alpha^(2*nu)*gamma(nu+0.5));
        data_matrix = 1/n*randn(n,n);

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
        G=diag(diag(G));
        A = Q*G*Q';
        figure(1)
        semilogy(diag(G))
        hold on 
        semilogy(diag(G)+mi*ones(n))
end

% lgA = log(diag(G+mi*eye(n)));
lgA = logm(A+mi*eye(n));

for N = 60000 : 10000 : 100000

    for j=1:N
        v = randn(n,1);
        esA_G(j) = v'*lgA*v;
        
        w = randsrc(n,1);
        esA_R(j) = w'*lgA*w;

        z = unifrnd(-1,1,n,1);
        z = sqrt(n)*z/norm(z);
        esA_U(j) = z'*lgA*z;
    end
    
    meanA_G = sum(esA_G(:))/N;
    varA_G(N/10000-5) = sum((esA_G-meanA_G).^2)/N;
    ctrlA_G = 2*norm(lgA,'fro')^2;

    meanA_R = sum(esA_R(:))/N;
    varA_R(N/10000-5) = sum((esA_R-meanA_R).^2)/N;
    ctrlA_R = 2*(norm(lgA,'fro')^2-norm(diag(lgA))^2);


    meanA_U = sum(esA_U(:))/N;
    varA_U(N/10000-5) = sum((esA_U-meanA_U).^2)/N;
    ctrlA_U = (2*n/(2+n))*(norm(lgA,'fro')^2-1/n*trace(lgA)^2);

end

figure(2)
plot(varA_G,'o-g')
hold on
plot(ctrlA_G*ones(1,5),'g')

figure(3)
plot(varA_R,'o-g')
hold on
plot(ctrlA_R*ones(1,5),'g')


figure(4)
plot(varA_U,'o-g')
hold on
plot(ctrlA_U*ones(1,5),'g')